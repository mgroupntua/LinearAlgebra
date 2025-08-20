using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;

using MGroup.Environments;

namespace MGroup.LinearAlgebra.Distributed.Overlapping
{
	/// <summary>
	/// All indexing data and functionality of <see cref="DistributedOverlappingIndexer"/>, but only for the local vector, 
	/// matrix, etc. that corresponds to a specific <see cref="ComputeNode"/>.
	/// </summary>
	internal class LocalIndexer
	{
		private Dictionary<int, int[]> commonEntriesWithNeighbors;

		internal LocalIndexer(ComputeNode node)
		{
			this.Node = node;
		}

		/// <summary>
		/// Neighboring <see cref="ComputeNode"/>s of this <see cref="Node"/> with local vectors that have at least 1 common 
		/// entry with the local vector of this <see cref="Node"/>.
		/// </summary>
		internal SortedSet<int> ActiveNeighborsOfNode { get; private set; }

		internal double[] InverseMultiplicities { get; private set; }

		internal ComputeNode Node { get; }

		internal int NumEntries { get; private set; }

		internal (int local, int remote) CountCommonEntries()
		{
			int local = 0;
			int remote = 0;
			foreach (var pair in commonEntriesWithNeighbors)
			{
				int neighborID = pair.Key;
				int[] commonEntries = pair.Value;
				if (this.Node.Cluster.Nodes.ContainsKey(neighborID))
				{
					local += commonEntries.Length;
				}
				else
				{
					remote += commonEntries.Length;
				}
			}
			return (local, remote);
		}

		//TODO: cache a buffer for sending and a buffer for receiving inside Indexer (lazily or not) and just return them. 
		//      Also provide an option to request newly initialized buffers. It may be better to have dedicated Buffer classes to
		//      handle all that logic (e.g. keeping allocated buffers in a LinkedList, giving them out & locking them, 
		//      freeing them in clients, etc.
		internal ConcurrentDictionary<int, double[]> CreateBuffersForAllToAllWithNeighbors()
		{
			//TODOMPI: dictionaries that contain per node values should be requested from the environment, which knows their
			//      type (Dictionary/ConcurrentDictionary), capacity and concurrency level.
			var buffers = new ConcurrentDictionary<int, double[]>();
			foreach (int neighborID in ActiveNeighborsOfNode)
			{
				buffers[neighborID] = new double[commonEntriesWithNeighbors[neighborID].Length];
			}
			return buffers;
		}

		internal LocalIndexer DeepCopy()
		{
			var clone = new LocalIndexer(this.Node);
			clone.NumEntries = this.NumEntries;
			clone.ActiveNeighborsOfNode = new SortedSet<int>(this.ActiveNeighborsOfNode);

			clone.InverseMultiplicities = new double[this.InverseMultiplicities.Length];
			Array.Copy(this.InverseMultiplicities, clone.InverseMultiplicities, this.InverseMultiplicities.Length);

			clone.commonEntriesWithNeighbors = new Dictionary<int, int[]>();
			foreach ((int nodeID, int[] data) in this.commonEntriesWithNeighbors)
			{
				var clonedData = new int[data.Length];
				Array.Copy(data, clonedData, data.Length);
				clone.commonEntriesWithNeighbors[nodeID] = clonedData;
			}

			return clone;
		}


		internal int[] GetCommonEntriesWithNeighbor(int neighborID) => commonEntriesWithNeighbors[neighborID];

		internal void Initialize(LocalIndexerDto indexingDto)
		{
			this.NumEntries = indexingDto.NumEntries;
			ActiveNeighborsOfNode = new SortedSet<int>(indexingDto.CommonEntriesOfNodeWithNeighbors.Keys);
			Debug.Assert(Node.Neighbors.IsSupersetOf(ActiveNeighborsOfNode));
			this.commonEntriesWithNeighbors = indexingDto.CommonEntriesOfNodeWithNeighbors;
			FindMultiplicities();
		}

		/// <summary>
		/// Copy data shallowly from <paramref name="other"/>.
		/// </summary>
		/// <param name="other"></param>
		internal void InitializeFrom(LocalIndexer other)
		{
			this.NumEntries = other.NumEntries;
			this.commonEntriesWithNeighbors = other.commonEntriesWithNeighbors;
			this.ActiveNeighborsOfNode = other.ActiveNeighborsOfNode;
			this.InverseMultiplicities = other.InverseMultiplicities;
		}

		private void FindMultiplicities()
		{
			var multiplicities = new int[NumEntries];
			for (int i = 0; i < NumEntries; ++i) multiplicities[i] = 1;
			foreach (int[] commonEntries in commonEntriesWithNeighbors.Values)
			{
				foreach (int i in commonEntries) multiplicities[i] += 1;
			}

			InverseMultiplicities = new double[NumEntries];
			for (int i = 0; i < NumEntries; ++i) InverseMultiplicities[i] = 1.0 / multiplicities[i];
		}
	}
}
