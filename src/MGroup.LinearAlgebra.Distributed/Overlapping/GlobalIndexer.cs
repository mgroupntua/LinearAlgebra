using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Runtime.Serialization;
using System.Text;

using MGroup.Environments;
using MGroup.MSolve.Discretization.Entities;

namespace MGroup.LinearAlgebra.Distributed.Overlapping
{
	/// <summary>
	/// In the distributed overlapping vector (and matrix) strategy, each entry of the global vector entry may appear in one or 
	/// more local vectors. This class is responsible for converting global indices to local ones.
	/// </summary>
	public class GlobalIndexer
    {
		/// <summary>
		/// Each entry of the List corresponds to one global index. Each global vector entry may appear in one or more local
		/// vectors. Thus for each global index we store a Dictionary&lt;int, int&gt;, where a) each key is the id of the
		/// ComputeNode corresponding to the local vector, b) each value is the index of the entry in the local vector.
		/// </summary>
		private readonly List<Dictionary<int, int>> globalToLocal;

		/// <summary>
		/// Each key of the dictionary is the index of the ComputeNode corresponding to a local vector. Each value is an array
		/// that stores the global indices that correspond to the local indices.
		/// </summary>
		private readonly Dictionary<int, int[]> localToGlobal;

		public GlobalIndexer(DistributedOverlappingIndexer mainIndexer)
		{
			Dictionary<int, DistributedOverlappingIndexer.Local> localIndexers = AllGatherLocalIndexers(mainIndexer);
			int numNodes = localIndexers.Count;

			// Allocate memory for global-to-local maps. 
			NumGlobalIndices = mainIndexer.CountUniqueEntries();
			globalToLocal = new List<Dictionary<int, int>>(NumGlobalIndices);
			for (int i = 0; i < NumGlobalIndices; i++)
			{
				globalToLocal.Add(new Dictionary<int, int>());
			}

			// Allocate memory for local-to-global maps and initialize to -1. 
			localToGlobal = new Dictionary<int, int[]>();
			for (int nodeID = 0; nodeID < numNodes; nodeID++)
			{
				DistributedOverlappingIndexer.Local localIndexer = localIndexers[nodeID];
				localToGlobal[nodeID] = new int[localIndexer.NumEntries];
				Array.Fill(localToGlobal[nodeID], -1);
			}

			// Iterate over the local vectors
			int globalIdx = -1;
			for (int nodeID = 0; nodeID < numNodes; nodeID++)
			{
				//ComputeNode node = mainIndexer.Environment.GetComputeNode(nodeID);
				DistributedOverlappingIndexer.Local localIndexer = localIndexers[nodeID];
				int numLocalEntries = localIndexers[nodeID].NumEntries;
				for (int localIdx = 0; localIdx < numLocalEntries; localIdx++)
				{
					if (localToGlobal[nodeID][localIdx] == -1) // This is an entry we have not met when processing another local vector
					{
						globalIdx++;
						globalToLocal[globalIdx].Add(nodeID, localIdx);
						localToGlobal[nodeID][localIdx] = globalIdx;

						// Do the same for other local vectors that contain this entry
						List<(int neighborNodeID, int neighborLocalIdx)> commonEntries =
							FindCommonEntriesInNeighbors(nodeID, localIdx, localIndexers);
						foreach ((int neighborNodeID, int neighborLocalIdx) in commonEntries)
						{
							globalToLocal[globalIdx].Add(neighborNodeID, neighborLocalIdx);
							localToGlobal[neighborNodeID][neighborLocalIdx] = globalIdx;
						}
					}
				}
			}

			// Check the mappings
			for (int i = 0; i < NumGlobalIndices; i++)
			{
				if (globalToLocal[i].Count < 1)
				{
					throw new Exception("The mapping from global-to-local indices could not be constructed.");
				}
			}
			for (int nodeID = 0; nodeID < numNodes; nodeID++)
			{
				DistributedOverlappingIndexer.Local localIndexer = localIndexers[nodeID];
				int[] map = localToGlobal[nodeID];
				for (int localIdx = 0; localIdx < map.Length; localIdx++)
				{
					if (map[localIdx] == -1)
					{
						throw new Exception("The mapping from local-to-global indices could not be constructed.");
					}
				}
			}
		}

		public int NumGlobalIndices { get; }

		public int FindGlobalIndexOf(int nodeID, int localIdx) => localToGlobal[nodeID][localIdx];

		/// <summary>
		/// Returns -1 if no such index exists
		/// </summary>
		/// <param name="globalIdx"></param>
		/// <param name="nodeID"></param>
		/// <returns></returns>
		/// <exception cref="ArgumentException">Invalid global index</exception>
		public int FindLocalIndexOf(int globalIdx, int nodeID)
		{
			if (globalIdx < 0 || globalIdx > globalToLocal.Count)
			{
				throw new ArgumentException(
					$"The are {globalToLocal.Count} global indices, but {globalIdx} was requested.");
			}

			Dictionary<int, int> nodeToLocal = globalToLocal[globalIdx];
			if (nodeToLocal.TryGetValue(nodeID, out int localIdx))
			{
				return localIdx;
			}
			else
			{
				return -1;
			}
		}

		public List<(int nodeID, int localIdx)> FindLocalIndicesOf(int globalIdx)
		{
			if (globalIdx < 0 || globalIdx > globalToLocal.Count)
			{
				throw new ArgumentException(
					$"The are {globalToLocal.Count} global indices, but {globalIdx} was requested.");
			}

			Dictionary<int, int> nodeToLocal = globalToLocal[globalIdx];
			var result = new List<(int nodeID, int localIdx)>(nodeToLocal.Count);
			foreach (var nodeId_localIdxPair in nodeToLocal)
			{
				result.Add((nodeId_localIdxPair.Key, nodeId_localIdxPair.Value));
			}

			return result;
		}

		private Dictionary<int, DistributedOverlappingIndexer.Local> AllGatherLocalIndexers(DistributedOverlappingIndexer mainIndexer)
		{
			Dictionary<int, LocalIndexerDto> transferedDtos = mainIndexer.Environment.AllGather(
				nodeID => new LocalIndexerDto(mainIndexer.GetLocalComponent(nodeID)));

			var localIndexers = new Dictionary<int, DistributedOverlappingIndexer.Local>();
			foreach (var nodeID_indexerDtoPair in transferedDtos)
			{
				localIndexers[nodeID_indexerDtoPair.Key] = nodeID_indexerDtoPair.Value.ToLocalIndexer(mainIndexer.Environment);
			}

			return localIndexers;
		}

		private List<(int neighborNodeID, int neighborLocalIdx)> FindCommonEntriesInNeighbors(int targetNodeID, 
			int targetLocalIdx, Dictionary<int, DistributedOverlappingIndexer.Local> localIndexers)
		{
			var result = new List<(int, int)>();
			DistributedOverlappingIndexer.Local targetLocalIndexer = localIndexers[targetNodeID];
			foreach (int neighborNodeID in targetLocalIndexer.ActiveNeighborsOfNode)
			{
				// Check if the target local index corresponds to a common entry with this neighbor
				int[] localIndices = targetLocalIndexer.GetCommonEntriesWithNeighbor(neighborNodeID);
				int pos = -1;
				for (int i = 0; i < localIndices.Length; i++)
				{
					if (localIndices[i] == targetLocalIdx)
					{
						pos = i;
						break;
					}
				}

				if (pos != -1)
				{
					int[] neighborLocalIndices = localIndexers[neighborNodeID].GetCommonEntriesWithNeighbor(targetNodeID);
					Debug.Assert(localIndices.Length == neighborLocalIndices.Length);
					result.Add((neighborNodeID, neighborLocalIndices[pos]));
				}
			}

			return result;
		}

		[Serializable]
		private class LocalIndexerDto
		{
			public LocalIndexerDto(DistributedOverlappingIndexer.Local localIndexer)
			{
				NodeID = localIndexer.Node.ID;
				NumEntries = localIndexer.NumEntries;
				CommonEntriesWithNeighbors = new Dictionary<int, int[]>();
				foreach (int neighborID in localIndexer.ActiveNeighborsOfNode)
				{
					CommonEntriesWithNeighbors[neighborID] = localIndexer.GetCommonEntriesWithNeighbor(neighborID);
				}
			}

			public Dictionary<int, int[]> CommonEntriesWithNeighbors { get; set; }

			public int NodeID { get; set; }

			public int NumEntries { get; set; }

			public DistributedOverlappingIndexer.Local ToLocalIndexer(IComputeEnvironment environment)
			{
				var localIndexer = new DistributedOverlappingIndexer.Local(environment.GetComputeNode(NodeID));
				localIndexer.Initialize(NumEntries, CommonEntriesWithNeighbors);
				return localIndexer;
			}
		}
	}
}
