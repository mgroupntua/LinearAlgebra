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
		private readonly Dictionary<int, LocalIndexer> localIndexers;
		private readonly object myLock = new();

		public DistributedOverlappingIndexer(IComputeEnvironment environment)
		{
			localIndexers = environment.CalcNodeData(
				n => new LocalIndexer(environment.GetComputeNode(n)));
			Environment = environment;
		}

		public IComputeEnvironment Environment { get; }

		public int NumGlobalIndices { get; private set; }

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

		public SortedSet<int> GetActiveNeighborIDs(int nodeID) => localIndexers[nodeID].ActiveNeighborsOfNode;

		public int[] GetCommonEntriesOfNodeWithNeighbor(int nodeID, int neighborID) 
			=> localIndexers[nodeID].GetCommonEntriesWithNeighbor(neighborID);

		public double[] GetInverseMultiplicities(int nodeID) => localIndexers[nodeID].InverseMultiplicities;

		public int GetNumLocalIndices(int nodeID) => localIndexers[nodeID].NumEntries;

		public void Initialize(Func<int, LocalIndexerDto> getLocalIndexingData)
		{
			Environment.DoPerNode(node =>
			{
				LocalIndexerDto localIndexingData = getLocalIndexingData(node);
				localIndexers[node].Initialize(localIndexingData);
			});
			CountUniqueEntries();
		}

		public bool IsCompatibleWith(IDistributedIndexer other) => this == other;

		public DistributedOverlappingIndexer ReplaceWithNewIndexer(Func<int, LocalIndexerDto> getLocalIndexingData)
		{
			var result = new DistributedOverlappingIndexer(Environment);
			Environment.DoPerNode(node =>
			{
				LocalIndexerDto localIndexingData = getLocalIndexingData(node);
				if (localIndexingData.Modified)
				{
					result.localIndexers[node].InitializeFrom(this.localIndexers[node]);
				}
				else
				{
					result.localIndexers[node].Initialize(localIndexingData);
				}
			});
			return result;
		}

		internal Dictionary<int, LocalIndexer> AllGatherLocalIndexers()
		{
			Dictionary<int, LocalIndexerDto> transferedDtos = Environment.AllGather(
				nodeID => new LocalIndexerDto(this.localIndexers[nodeID])
			);

			var result = new Dictionary<int, LocalIndexer>();
			foreach (var nodeID_indexerDtoPair in transferedDtos)
			{
				result[nodeID_indexerDtoPair.Key] = nodeID_indexerDtoPair.Value.ToLocalIndexer(Environment);
			}

			return result;
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
