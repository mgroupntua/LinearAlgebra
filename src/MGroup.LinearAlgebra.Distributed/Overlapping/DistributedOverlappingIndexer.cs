using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;
using MGroup.Environments;
using MGroup.LinearAlgebra.Exceptions;
using MGroup.LinearAlgebra.Matrices;
using MGroup.LinearAlgebra.Vectors;

namespace MGroup.LinearAlgebra.Distributed.Overlapping
{
	/// <summary>
	/// Manages the indices for a <see cref="DistributedOverlappingVector"/>, <see cref="DistributedOverlappingMatrix{TMatrix}"/>, 
	/// etc. Supports multiple local vectors, each of which may have none, some or all its entries in common with other local 
	/// vectors/matrices. Specifies the relationships between these common entries. When dealing with multiple distributed  
	/// vectors that have the same indexing pattern, reuse the same instance of <see cref="DistributedOverlappingIndexer"/>.
	/// </summary>
	/// <remarks>
	/// In interface problems of PSM and related DDMs, local vectors have all their entries in common with other local vectors, 
	/// since only boundary dofs take part in the interface problem. In GSI or a GSI-like treatment of other DDMs' coarse 
	/// problems, local vectors have some of their entries (boundary dofs) in common with other local vectors, while the rest
	/// entries (internal dofs) are unique for each local vector. 
	/// </remarks>
	public class DistributedOverlappingIndexer : IDistributedIndexer
	{
		private readonly object myLock = new();
		private Dictionary<int, LocalIndexer> localIndexers;
		private GlobalIndexer globalIndexer;

		public DistributedOverlappingIndexer(IComputeEnvironment environment)
		{
			Environment = environment;
		}

		public IComputeEnvironment Environment { get; }

		public int NumGlobalIndices { get; private set; }

		public void CheckGlobalIndex1D(int index)
		{
			if (index < 0 || index >= NumGlobalIndices)
			{
				throw new IndexOutOfRangeException($"The index must be in the range [0, {NumGlobalIndices}), but was {index}");
			}
		}

		internal void CheckGlobalIndex2D(int rowIdx, int colIdx)
		{
			if (rowIdx < 0 || rowIdx >= NumGlobalIndices)
			{
				throw new IndexOutOfRangeException(
					$"The row index must be in the range [0, {NumGlobalIndices}), but was {rowIdx}");
			}

			if (colIdx < 0 || colIdx >= NumGlobalIndices)
			{
				throw new IndexOutOfRangeException(
					$"The column index must be in the range [0, {NumGlobalIndices}), but was {colIdx}");
			}
		}

		/// <summary>
		/// Counts how many vector (or matrix) entries are common with other compute nodes. 
		/// </summary>
		/// <returns>
		/// "local" = how many entries in common with vectors (or matrices) in the same cluster.
		/// "remote" = how many entries in common with vectors (or matrices) in other clusters.
		/// </returns>
		public (int local, int remote) CountCommonEntriesOfNodeWithNeighbors(int nodeID) 
			=> localIndexers[nodeID].CountCommonEntries();

		public ConcurrentDictionary<int, double[]> CreateBuffersForAllToAllWithNeighbors(int nodeID) 
			=> localIndexers[nodeID].CreateBuffersForAllToAllWithNeighbors();

		public DistributedOverlappingIndexer DeepCopy()
		{
			var clone = new DistributedOverlappingIndexer(this.Environment);
			Environment.DoPerNode(node =>
			{
				clone.localIndexers[node] = this.localIndexers[node].DeepCopy();
			});
			clone.NumGlobalIndices = this.NumGlobalIndices;
			return clone;
		}

		public int FindGlobalIndexOf(int nodeID, int localIdx) 
			=> CreateGlobalIndexerIfMissing().FindGlobalIndexOf(nodeID, localIdx);

		/// <summary>
		/// Returns the global index corresponding to a local index or -1 if no such entry exists.
		/// </summary>
		/// <param name="globalIdx">The global index of the entry.</param>
		/// <param name="nodeID">The id of the local vector/node</param>
		/// <returns>See summary</returns>
		/// <exception cref="ArgumentException">Invalid global index</exception>
		public int FindLocalIndexOf(int globalIdx, int nodeID)
			=> CreateGlobalIndexerIfMissing().FindLocalIndexOf(globalIdx, nodeID);

		/// <summary>
		/// Returns a dictionary where: a) keys are the ids of the nodes (1 node -> 1 local vector) containing
		/// <paramref name="globalIdx"/>, b) values are the local indices for the corresponding nodes.
		/// </summary>
		public IReadOnlyDictionary<int, int> FindLocalIndicesOf(int globalIdx) 
			=> CreateGlobalIndexerIfMissing().FindLocalIndicesOf(globalIdx);

		public SortedSet<int> GetActiveNeighborIDs(int nodeID) => localIndexers[nodeID].ActiveNeighborsOfNode;

		public int[] GetCommonEntriesOfNodeWithNeighbor(int nodeID, int neighborID) 
			=> localIndexers[nodeID].GetCommonEntriesWithNeighbor(neighborID);

		public double[] GetInverseMultiplicities(int nodeID) => localIndexers[nodeID].InverseMultiplicities;

		public int GetNumLocalIndices(int nodeID) => localIndexers[nodeID].NumIndices;

		public void Initialize(Func<int, LocalIndexerDto> getLocalIndexingData)
		{
			localIndexers = Environment.CalcNodeData(nodeID =>
			{
				LocalIndexerDto dto = getLocalIndexingData(nodeID);
				return new LocalIndexer(Environment.GetComputeNode(nodeID), dto.CommonEntriesOfNodeWithNeighbors, dto.NumIndices);
			});
			CountUniqueEntries();
		}

		public bool IsCompatibleWith(IDistributedIndexer other) => this == other;


		public DistributedOverlappingIndexer ReuseAsBasisForNewIndexer(Func<int, LocalIndexerDto> getLocalIndexingData)
		{
			var result = new DistributedOverlappingIndexer(Environment);
			result.localIndexers = Environment.CalcNodeData(nodeID =>
			{
				LocalIndexerDto dto = getLocalIndexingData(nodeID);
				if (dto.Modified)
				{
					return new LocalIndexer(
						Environment.GetComputeNode(nodeID), dto.CommonEntriesOfNodeWithNeighbors, dto.NumIndices);
				}
				else
				{
					return this.localIndexers[nodeID];
				}
			});

			result.CountUniqueEntries();
			return result;
		}

		internal Dictionary<int, LocalIndexer> AllGatherLocalIndexers()
		{
			Dictionary<int, LocalIndexerDto> transferedDtos = Environment.AllGather(
				nodeID => LocalIndexerDto.CreateForSerialization(this.localIndexers[nodeID])
			);

			var result = new Dictionary<int, LocalIndexer>();
			foreach (var nodeID_indexerDtoPair in transferedDtos)
			{
				result[nodeID_indexerDtoPair.Key] = nodeID_indexerDtoPair.Value.ToLocalIndexer(Environment);
			}

			return result;
		}

		private GlobalIndexer CreateGlobalIndexerIfMissing()
		{
			if (globalIndexer == null)
			{
				lock (myLock)
				{
					if (globalIndexer == null) // in case another thread created it before this thread got the lock
					{
						globalIndexer = new GlobalIndexer(localIndexers, NumGlobalIndices);
					}
				}
			}

			return globalIndexer;
		}

		private void CountUniqueEntries()
		{
			Dictionary<int, double> countPerNode = Environment.CalcNodeData(node =>
			{
				double[] inverseMultiplicities = localIndexers[node].InverseMultiplicities;

				double localCount = 0.0;
				for (int i = 0; i < inverseMultiplicities.Length; ++i)
				{
					localCount += inverseMultiplicities[i];
				}

				return localCount;
			});
			double globalCount = Environment.AllReduceSum(countPerNode);
			NumGlobalIndices = (int)Math.Round(globalCount);
		}
	}
}
