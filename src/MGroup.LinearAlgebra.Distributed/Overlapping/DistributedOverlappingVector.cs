using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;
using MGroup.LinearAlgebra.Vectors;
using MGroup.Environments;
using MGroup.LinearAlgebra.Distributed.LinearAlgebraExtensions;
using System.Collections.Concurrent;
using MGroup.MSolve.Solution.LinearSystem;

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
	public class DistributedOverlappingVector : IGlobalVector
	{
		private ConcurrentDictionary<int, (ConcurrentDictionary<int, double[]> send, ConcurrentDictionary<int, double[]> recv)>	cachedBuffers = 
			new ConcurrentDictionary<int, (ConcurrentDictionary<int, double[]> send, ConcurrentDictionary<int, double[]> recv)>();

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

		public IDictionary<int, Vector> LocalVectors { get; }

		public bool CheckForCompatibility { get; set; } = true;

		public void AxpyIntoThis(IGlobalVector otherVector, double otherCoefficient)
		{
			DistributedOverlappingVector otherDistributed = Indexer.CheckCompatibleVector(otherVector);
			AxpyIntoThis(otherDistributed, otherCoefficient);
		}

		public void AxpyIntoThis(DistributedOverlappingVector otherVector, double otherCoefficient)
		{
			Indexer.CheckCompatibleVector(otherVector);
			Environment.DoPerNode(
				node => this.LocalVectors[node].AxpyIntoThis(otherVector.LocalVectors[node], otherCoefficient)
			);
		}

		public void Clear()
		{
			Environment.DoPerNode(node => LocalVectors[node].Clear());
		}

		IGlobalVector IGlobalVector.Copy() => Copy(); //TODO: Copy can be expressed with CreateZero() and CopyFrom().

		public DistributedOverlappingVector Copy()
		{
			Dictionary<int, Vector> localVectorsCloned =
				Environment.CalcNodeData(node => LocalVectors[node].Copy());
			return new DistributedOverlappingVector(Indexer, localVectorsCloned);
		}

		public void CopyFrom(IGlobalVector otherVector)
		{
			DistributedOverlappingVector otherDistributed = Indexer.CheckCompatibleVector(otherVector);
			CopyFrom(otherDistributed);
		}

		public void CopyFrom(DistributedOverlappingVector otherVector)
		{
			Indexer.CheckCompatibleVector(otherVector);
			Environment.DoPerNode(node => this.LocalVectors[node].CopyFrom(otherVector.LocalVectors[node]));
		}

		IGlobalVector IGlobalVector.CreateZero() => CreateZero();

		public DistributedOverlappingVector CreateZero()
		{
			var result = new DistributedOverlappingVector(Indexer);
			result.CacheSendRecvBuffers = this.CacheSendRecvBuffers;
			return result;
		}

		public double DotProduct(IGlobalVector otherVector)
		{
			DistributedOverlappingVector otherDistributed = Indexer.CheckCompatibleVector(otherVector);
			return DotProduct(otherDistributed);
		}

		/// <summary>
		/// See <see cref="IGlobalVector.DotProduct(IGlobalVector)"/>.
		/// </summary>
		/// <remarks>
		/// Warning: This does not work correctly if 2 local vectors have different values at the same common entry. In such 
		/// cases make, perhaps <see cref="SumOverlappingEntries"/> may be of use.
		/// </remarks>
		public double DotProduct(DistributedOverlappingVector otherVector)
		{
			Indexer.CheckCompatibleVector(otherVector);
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

		public int Length()
		{
			Dictionary<int, double> localLengths = Environment.CalcNodeData(node =>
			{
				double[] inverseMultiplicities = Indexer.GetLocalComponent(node).InverseMultiplicities;

				double localLegth = 0.0;
				for (int i = 0; i < inverseMultiplicities.Length; ++i)
				{
					localLegth += inverseMultiplicities[i];
				}

				return localLegth;
			});

			return (int)Math.Round(Environment.AllReduceSum(localLengths));
		}

		public void LinearCombinationIntoThis(
			double thisCoefficient, IGlobalVector otherVector, double otherCoefficient)
		{
			DistributedOverlappingVector otherDistributed = Indexer.CheckCompatibleVector(otherVector);
			LinearCombinationIntoThis(thisCoefficient, otherDistributed, otherCoefficient);
		}

		public void LinearCombinationIntoThis(
			double thisCoefficient, DistributedOverlappingVector otherVector, double otherCoefficient)
		{
			Indexer.CheckCompatibleVector(otherVector);
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

		public void ScaleIntoThis(double scalar)
		{
			Environment.DoPerNode(node => LocalVectors[node].ScaleIntoThis(scalar));
		}

		public void SetAll(double value)
		{
			Environment.DoPerNode(node => LocalVectors[node].SetAll(value));
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

		private (ConcurrentDictionary<int, double[]> sendValues, ConcurrentDictionary<int, double[]> recvValues) 
			GetSendRecvBuffers(int nodeID)
		{
			#region debug
			//CacheSendRecvBuffers = false;
			#endregion
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
