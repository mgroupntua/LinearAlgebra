using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;
using MGroup.LinearAlgebra.Vectors;
using MGroup.Environments;
using System.Collections.Concurrent;
using MGroup.LinearAlgebra.Commons;
using MGroup.LinearAlgebra.Distributed.Utilities;
using MGroup.LinearAlgebra.Exceptions;
using MGroup.LinearAlgebra.Reduction;
using DotNumerics.FortranLibrary;

using static MGroup.LinearAlgebra.Distributed.Overlapping.CompatibilityUtilities;

//TODOMPI: this class will be mainly used for iterative methods. Taking that into account, make optimizations. E.g. work arrays
//      used as buffers for MPI communication can be reused across vectors, instead of each vector allocating/freeing identical 
//      buffers. Such functionality can be included in the indexer, which is shared across vectors/matrices.
//TODOMPI: should this class have a Length property? It seems important for many linear algebra dimension matching checks, but 
//      it will probably require significant communication. Furthermore, these checks can probably depend on polymorphic methods
//      exposed by the vectors & matrix classes, which will check matching dimensions between matrix-vector or vector-vector.
//      E.g. Adding 2 vectors requires that they have the same length. Vector will check exactly that, and possibly expose a 
//      PatternMatchesForLinearCombo(other) method. DistributedVector however will check that they have the same indexers, 
//      without any need to communicate, only to find the total length. If I do provide such a property, it should be accessed 
//      from the indexer (which must be 1 object for all compute nodes). The indexer should lazily calculate it, store it
//      internally and update it whenever the connectivity changes. Or just prohibit changing the connectivity. Calculating it
//      will be similar to the dot product: sum the number of internal and boundary entries in each local node (divide the 
//      boundary entries over the multiplicities resulting in fractional number), reduce the double result from node and finally
//      round it to the nearest integer (and pray the precision errors are negligible).
namespace MGroup.LinearAlgebra.Distributed.Overlapping
{
	public class DistributedOverlappingVector : IVector
	{
		private ConcurrentDictionary<int, (ConcurrentDictionary<int, double[]> send, ConcurrentDictionary<int, double[]> recv)>	cachedBuffers = 
			new ConcurrentDictionary<int, (ConcurrentDictionary<int, double[]> send, ConcurrentDictionary<int, double[]> recv)>();

		private GlobalIndexer globalIndexer;

		public DistributedOverlappingVector(DistributedOverlappingIndexer indexer)
		{
			this.Indexer = indexer;
			this.Environment = indexer.Environment;
			this.LocalVectors = Environment.CalcNodeData(
				node => Vector.CreateZero(indexer.GetLocalComponent(node).NumEntries));
		}

		public DistributedOverlappingVector(DistributedOverlappingIndexer indexer, IDictionary<int, Vector> localVectors)
		{
			this.Indexer = indexer;
			this.Environment = indexer.Environment;
			this.LocalVectors = localVectors;
		}

		public DistributedOverlappingVector(DistributedOverlappingIndexer indexer, Func<int, Vector> createLocalVector)
		{
			this.Indexer = indexer;
			this.Environment = indexer.Environment;
			this.LocalVectors = Environment.CalcNodeData(createLocalVector);
		}

		public bool CacheSendRecvBuffers { get; set; } = false;

		public IComputeEnvironment Environment { get; }

		public DistributedOverlappingIndexer Indexer { get; }

		public int Length => Indexer.NumUniqueEntries;

		public IDictionary<int, Vector> LocalVectors { get; }

		public double this[int index]
		{
			get 
			{
				IReadOnlyDictionary<int, int> localIndices = FindLocalIndicesFromGlobal(index);
				foreach ((int nodeID, int localIdx) in localIndices)
				{
					if (LocalVectors.TryGetValue(nodeID, out Vector localVector))
					{
						return localVector[localIdx]; // Return the first to be found.
					}
				}

				throw new Exception("This should not have happened. The distributed vector is not created correctly.");
			}
		}

		public void AddToIndex(int index, double value)
		{
			IReadOnlyDictionary<int, int> localIndices = FindLocalIndicesFromGlobal(index);
			foreach ((int nodeID, int localIdx) in localIndices)
			{
				if (LocalVectors.TryGetValue(nodeID, out Vector localVector))
				{
					localVector.AddToIndex(localIdx, value); // Set all instances
				}
			}
		}

		public bool AreOverlappingEntriesEqual(double tolerance)
		{
			var comparer = new ValueComparer(tolerance);

			Dictionary<int, AllToAllNodeData<double>> dataPerNode = ExchangeOverlappingEntries();

			// Add the common entries of neighbors back to the original local vector.
			Func<int, bool> checkLocalSubvector = nodeID =>
			{
				ComputeNode node = Environment.GetComputeNode(nodeID);
				Vector localVector = LocalVectors[nodeID];
				DistributedOverlappingIndexer.Local localIndexer = Indexer.GetLocalComponent(nodeID);

				IDictionary<int, double[]> recvValues = dataPerNode[nodeID].recvValues;
				foreach (int neighborID in localIndexer.ActiveNeighborsOfNode)
				{
					int[] commonEntries = localIndexer.GetCommonEntriesWithNeighbor(neighborID);
					Vector localValues = localVector.GetSubvector(commonEntries);
					double[] neighborValues = recvValues[neighborID];
					Debug.Assert(localValues.Length == neighborValues.Length);
					if (node.ID < neighborID) // Make sure the comparisons are done with identical arguments for both compute nodes
					{
						for (int i = 0; i < localValues.Length; ++i)
						{
							if (!comparer.AreEqual(localValues[i], neighborValues[i]))
							{
								return false;
							}
						}
					}
					else
					{
						for (int i = 0; i < localValues.Length; ++i)
						{
							if (!comparer.AreEqual(neighborValues[i], localValues[i]))
							{
								return false;
							}
						}
					}
				}
				return true;
			};
			Dictionary<int, bool> localResults = Environment.CalcNodeData(checkLocalSubvector);
			return Environment.AllReduceAnd(localResults);
		}

		public IVector Axpy(IVectorView otherVector, double otherCoefficient)
		{
			DistributedOverlappingVector result = Copy();
			result.AxpyIntoThis(otherVector, otherCoefficient);
			return result;
		}

		public void AxpyIntoThis(IVectorView otherVector, double otherCoefficient) 
			=> AxpyIntoThis(CastToDistributed(otherVector), otherCoefficient);

		public void AxpyIntoThis(DistributedOverlappingVector otherVector, double otherCoefficient)
		{
			CheckSameFormat(this, otherVector);
			Environment.DoPerNode(
				node => this.LocalVectors[node].AxpyIntoThis(otherVector.LocalVectors[node], otherCoefficient)
			);
		}

		public void Clear()
		{
			Environment.DoPerNode(node => LocalVectors[node].Clear());
		}

		IVector IVectorView.Copy(bool copyIndexingData = false) => Copy(copyIndexingData); //TODO: Copy can be expressed with CreateZero() and CopyFrom().

		public DistributedOverlappingVector Copy(bool copyIndexingData = false)
		{
			var indexerCloned = copyIndexingData ? Indexer.DeepCopy() : Indexer;
			Dictionary<int, Vector> localVectorsCloned =
				Environment.CalcNodeData(node => LocalVectors[node].Copy());
			return new DistributedOverlappingVector(indexerCloned, localVectorsCloned);
		}

		public void CopyFrom(IVectorView otherVector) => CopyFrom(CastToDistributed(otherVector));

		public void CopyFrom(DistributedOverlappingVector otherVector)
		{
			CheckSameFormat(this, otherVector);
			Environment.DoPerNode(node => this.LocalVectors[node].CopyFrom(otherVector.LocalVectors[node]));
		}

		public double[] CopyToArray()
		{
			CreateGlobalIndexerIfMissing();
			var result = new double[Length];
			Environment.DoPerNodeSerially(nodeID =>
			{
				Vector localVector = LocalVectors[nodeID];
				for (int localIdx = 0; localIdx < localVector.Length; localIdx++)
				{
					int globalIdx = globalIndexer.FindGlobalIndexOf(nodeID, localIdx);
					result[globalIdx] = localVector[localIdx];
				}
			});
			return result;
		}

		IVector IVectorView.CreateZeroVectorWithSameFormat() => CreateZeroVectorWithSameFormat();

		public DistributedOverlappingVector CreateZeroVectorWithSameFormat()
		{
			var result = new DistributedOverlappingVector(Indexer);
			result.CacheSendRecvBuffers = this.CacheSendRecvBuffers;
			return result;
		}

		public IVector DoEntrywise(IVectorView other, Func<double, double, double> binaryOperation)
		{
			DistributedOverlappingVector result = Copy();
			result.DoEntrywiseIntoThis(other, binaryOperation);
			return result;
		}

		public void DoEntrywiseIntoThis(IVectorView vector, Func<double, double, double> binaryOperation)
			=> DoEntrywiseIntoThis(CastToDistributed(vector), binaryOperation);

		public void DoEntrywiseIntoThis(DistributedOverlappingVector other, Func<double, double, double> binaryOperation)
		{
			CheckSameFormat(this, other);
			Environment.DoPerNode(node => this.LocalVectors[node].DoEntrywiseIntoThis(other, binaryOperation));
		}

		public IVector DoToAllEntries(Func<double, double> unaryOperation)
		{
			DistributedOverlappingVector result = Copy();
			result.DoToAllEntriesIntoThis(unaryOperation);
			return result;
		}

		public void DoToAllEntriesIntoThis(Func<double, double> unaryOperation)
		{
			Environment.DoPerNode(node => this.LocalVectors[node].DoToAllEntriesIntoThis(unaryOperation));
		}

		public double DotProduct(IVectorView otherVector) => DotProduct(CastToDistributed(otherVector));

		/// <summary>
		/// See <see cref="IGlobalVector.DotProduct(IGlobalVector)"/>.
		/// </summary>
		/// <remarks>
		/// Warning: This does not work correctly if 2 local vectors have different values at the same common entry. In such 
		/// cases make, perhaps <see cref="SumOverlappingEntries"/> may be of use.
		/// </remarks>
		public double DotProduct(DistributedOverlappingVector otherVector)
		{
			CheckSameFormat(this, otherVector);
			Func<int, double> calcLocalDot = node =>
			{
				Vector thisLocalVector = this.LocalVectors[node];
				Vector otherLocalVector = otherVector.LocalVectors[node];
				double[] inverseMultiplicities = Indexer.GetLocalComponent(node).InverseMultiplicities;

				double dotLocal = 0.0;
				int length = thisLocalVector.Length;
				for (int i = 0; i < length; ++i)
				{
					dotLocal += thisLocalVector[i] * otherLocalVector[i] * inverseMultiplicities[i];
				}

				return dotLocal;
			};

			Dictionary<int, double> dotPerNode = Environment.CalcNodeData(calcLocalDot);
			return Environment.AllReduceSum(dotPerNode);
		}

		public bool Equals(IIndexable1D other, double tolerance = 1E-7)
		{
			if (other is DistributedOverlappingVector casted)
			{
				return this.Equals(casted, tolerance);
			}

			return false;
		}

		public bool Equals(DistributedOverlappingVector other, double tolerance = 1E-7)
		{
			if (!this.Indexer.IsCompatibleWith(other.Indexer))
			{
				return false;
			}

			Dictionary<int, bool> flags = Environment.CalcNodeData(
					node => this.LocalVectors[node].Equals(other.LocalVectors[node], tolerance));
			return Environment.AllReduceAnd(flags);
		}

		public bool HasSameFormat(IIndexable1D other)
		{
			if (other is DistributedOverlappingVector casted)
			{
				return this.Indexer.IsCompatibleWith(casted.Indexer);
			}

			return false;
		}

		public IVector LinearCombination(double thisCoefficient, IVectorView otherVector, double otherCoefficient)
		{
			DistributedOverlappingVector result = Copy();
			result.LinearCombinationIntoThis(thisCoefficient, otherVector, otherCoefficient);
			return result;
		}

		public void LinearCombinationIntoThis(
			double thisCoefficient, IVectorView otherVector, double otherCoefficient)
			=> LinearCombinationIntoThis(thisCoefficient, CastToDistributed(otherVector), otherCoefficient);

		public void LinearCombinationIntoThis(
			double thisCoefficient, DistributedOverlappingVector otherVector, double otherCoefficient)
		{
			CheckSameFormat(this, otherVector);
			Environment.DoPerNode(
				node => this.LocalVectors[node].LinearCombinationIntoThis(
					thisCoefficient, otherVector.LocalVectors[node], otherCoefficient)
			);
		}

		public double Norm2()
		{
			Func<int, double> calcLocalDot = node =>
			{
				Vector localVector = this.LocalVectors[node];
				double[] inverseMultiplicities = Indexer.GetLocalComponent(node).InverseMultiplicities;

				double dotLocal = 0.0;
				for (int i = 0; i < localVector.Length; ++i)
				{
					dotLocal += localVector[i] * localVector[i] * inverseMultiplicities[i];
				}

				return dotLocal;
			};

			Dictionary<int, double> dotPerNode = Environment.CalcNodeData(calcLocalDot);
			return Math.Sqrt(Environment.AllReduceSum(dotPerNode));
		}

		//TODOMPI: A ReduceOverlappingEntries(IReduction), which would cover sum and regularization would be more useful. 
		//      However the implementation should not be slower than the current SumOverlappingEntries(), since that is a very
		//      important operation.
		//TODOMPI: Test this
		/// <summary>
		/// Gathers the entries of remote vectors that correspond to the boundary entries of the local vectors and regularizes 
		/// them, meaning each of these entries is divided via the sum of corresponding entries over all local vectors. 
		/// Therefore, the resulting local vectors will not have the same values at their corresponding overlapping entries.
		/// </summary>
		/// <remarks>
		/// Requires communication between compute nodes:
		/// Each compute node sends its boundary entries to the neighbors that are assiciated with these entries. 
		/// Each neighbor receives only the entries it has in common.
		/// </remarks>
		public void RegularizeOverlappingEntries()
		{
			// Sum the values of overlapping entries in a different vector.
			DistributedOverlappingVector reducedVector = Copy();
			reducedVector.SumOverlappingEntries();

			// Divide the values of overlapping entries via their sums.
			Action<int> regularizeLocalVectors = nodeID =>
			{
				ComputeNode node = Environment.GetComputeNode(nodeID);
				DistributedOverlappingIndexer.Local localIndexer = Indexer.GetLocalComponent(nodeID);
				Vector orginalLocalVector = this.LocalVectors[nodeID];
				Vector reducedLocalVector = reducedVector.LocalVectors[nodeID];

				for (int i = 0; i < localIndexer.NumEntries; ++i)
				{
					//TODO: This assumes that all entries with multiplicity > 1 are overlapping and must be regularized. 
					//      Is that always a correct assumption?
					//TODO: Perhaps some tolerance should be used or the original int[] Multiplicities.
					if (localIndexer.InverseMultiplicities[i] < 1.0)
					{
						orginalLocalVector[i] /= reducedLocalVector[i];
					}
				}
			};
			Environment.DoPerNode(regularizeLocalVectors);
		}

		public IVector Scale(double scalar)
		{
			DistributedOverlappingVector result = Copy();
			result.ScaleIntoThis(scalar);
			return result;
		}

		public void ScaleIntoThis(double scalar)
		{
			Environment.DoPerNode(node => LocalVectors[node].ScaleIntoThis(scalar));
		}

		public void Set(int index, double value)
		{
			IReadOnlyDictionary<int, int> localIndices = FindLocalIndicesFromGlobal(index);
			foreach ((int nodeID, int localIdx) in localIndices)
			{
				if (LocalVectors.TryGetValue(nodeID, out Vector localVector))
				{
					localVector[localIdx] = value; // Set all instances
				}
			}
		}

		public void SetAll(double value)
		{
			Environment.DoPerNode(node => LocalVectors[node].SetAll(value));
		}

		/// <summary>
		/// Gathers the entries of remote vectors that correspond to the boundary entries of the local vectors and sums them.
		/// As a result, the overlapping entries of each local vector will have the same values. These values are the same
		/// as the ones we would have if a global vector was created by assembling the local vectors.
		/// </summary>
		/// <remarks>
		/// Requires communication between compute nodes:
		/// Each compute node sends its boundary entries to the neighbors that are assiciated with these entries. 
		/// Each neighbor receives only the entries it has in common.
		/// </remarks>
		public void SumOverlappingEntries()
		{
			Dictionary<int, AllToAllNodeData<double>> dataPerNode = ExchangeOverlappingEntries();

			// Add the common entries of neighbors back to the original local vector.
			Action<int> sumLocalSubvectors = nodeID =>
			{
				ComputeNode node = Environment.GetComputeNode(nodeID);
				Vector localVector = LocalVectors[nodeID];
				DistributedOverlappingIndexer.Local localIndexer = Indexer.GetLocalComponent(nodeID);

				IDictionary<int, double[]> recvValues = dataPerNode[nodeID].recvValues;
				foreach (int neighborID in localIndexer.ActiveNeighborsOfNode)
				{
					int[] commonEntries = localIndexer.GetCommonEntriesWithNeighbor(neighborID);
					var rv = Vector.CreateFromArray(recvValues[neighborID]);
					localVector.AddIntoThisNonContiguouslyFrom(commonEntries, rv);
				}
			};
			Environment.DoPerNode(sumLocalSubvectors);
		}

		private void CreateGlobalIndexerIfMissing()
		{
			if (globalIndexer == null)
			{
				lock (globalIndexer)
				{
					if (globalIndexer == null) // in case another thread created it before this thread got the lock
					{
						globalIndexer = new GlobalIndexer(Indexer);
					}
				}
			}
		}

		private IReadOnlyDictionary<int, int> FindLocalIndicesFromGlobal(int globalIdx)
		{
			CreateGlobalIndexerIfMissing();
			globalIndexer.CheckGlobalIndex1D(globalIdx);
			IReadOnlyDictionary<int, int> localIndices = globalIndexer.FindLocalIndicesOf(globalIdx);
			if (localIndices.Count == 0)
			{
				throw new Exception("This should not have happened. The distributed vector is not created correctly.");
			}

			return localIndices;
		}

		private Dictionary<int, AllToAllNodeData<double>> ExchangeOverlappingEntries()
		{
			// Prepare the boundary entries of each node before communicating them to its neighbors.
			Func<int, AllToAllNodeData<double>> prepareLocalData = nodeID =>
			{
				ComputeNode node = Environment.GetComputeNode(nodeID);
				Vector localVector = LocalVectors[nodeID];
				DistributedOverlappingIndexer.Local localIndexer = Indexer.GetLocalComponent(nodeID);

				// Find the common entries (to send and receive) of this node with each of its neighbors
				var transferData = new AllToAllNodeData<double>();
				(transferData.sendValues, transferData.recvValues) = GetSendRecvBuffers(nodeID);
				foreach (int neighborID in localIndexer.ActiveNeighborsOfNode)
				{
					int[] commonEntries = localIndexer.GetCommonEntriesWithNeighbor(neighborID);
					var sv = Vector.CreateFromArray(transferData.sendValues[neighborID]);
					sv.CopyNonContiguouslyFrom(localVector, commonEntries);
				}

				return transferData;
			};
			var dataPerNode = Environment.CalcNodeData(prepareLocalData);

			// Perform AllToAll to exchange the common boundary entries of each node with its neighbors.
			Environment.NeighborhoodAllToAll(dataPerNode, true);

			return dataPerNode;
		}

		private (ConcurrentDictionary<int, double[]> sendValues, ConcurrentDictionary<int, double[]> recvValues) 
			GetSendRecvBuffers(int nodeID)
		{
			if (CacheSendRecvBuffers)
			{
				bool isCached = cachedBuffers.TryGetValue(nodeID, 
					out (ConcurrentDictionary<int, double[]> send, ConcurrentDictionary<int, double[]> recv) buffers);
				if (!isCached)
				{
					DistributedOverlappingIndexer.Local localIndexer = Indexer.GetLocalComponent(nodeID);
					buffers = (localIndexer.CreateBuffersForAllToAllWithNeighbors(), 
						localIndexer.CreateBuffersForAllToAllWithNeighbors());
					cachedBuffers[nodeID] = buffers;
				}
				else
				{
					// No need to clear them as they will be overwritten.
					//foreach (double[] buffer in buffers.send.Values)
					//{
					//	Array.Clear(buffer, 0, buffer.Length);
					//}
					//foreach (double[] buffer in buffers.recv.Values)
					//{
					//	Array.Clear(buffer, 0, buffer.Length);
					//}
				}
				return buffers;
			}
			else
			{
				DistributedOverlappingIndexer.Local localIndexer = Indexer.GetLocalComponent(nodeID);
				ConcurrentDictionary<int, double[]> sendValues = localIndexer.CreateBuffersForAllToAllWithNeighbors();
				ConcurrentDictionary<int, double[]> recvValues = localIndexer.CreateBuffersForAllToAllWithNeighbors();
				return (sendValues, recvValues);
			}
		}
	}
}
