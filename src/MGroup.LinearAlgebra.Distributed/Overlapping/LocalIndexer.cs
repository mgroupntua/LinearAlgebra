using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;

using MGroup.Environments;

namespace MGroup.LinearAlgebra.Distributed.Overlapping
{
	/// <summary>
	/// All indexing data and functionality of <see cref="DistributedOverlappingIndexer"/>, but only for the local vector / 
	/// matrix / etc. that corresponds to a specific <see cref="ComputeNode"/>.
	/// </summary>
	internal class LocalIndexer
	{
		private readonly Dictionary<int, int[]> commonEntriesWithNeighbors;

		internal LocalIndexer(ComputeNode node, SortedSet<int> activeNeighbors, 
			Dictionary<int, int[]> commonEntriesWithNeighbors, double[] inverseMultiplicities, int numIndices)
		{
			this.Node = node;
			this.NumIndices = numIndices;
			this.commonEntriesWithNeighbors = commonEntriesWithNeighbors;
			this.ActiveNeighborsOfNode = activeNeighbors;
			this.InverseMultiplicities = inverseMultiplicities;
		}

		internal LocalIndexer(ComputeNode node, Dictionary<int, int[]> commonEntriesWithNeighbors, int numIndices)
		{
			this.Node = node;
			this.NumIndices = numIndices;
			this.commonEntriesWithNeighbors = commonEntriesWithNeighbors;
			ActiveNeighborsOfNode = new SortedSet<int>(commonEntriesWithNeighbors.Keys);
			Debug.Assert(Node.Neighbors.IsSupersetOf(ActiveNeighborsOfNode));
			this.InverseMultiplicities = FindMultiplicities();
		}

		/// <summary>
		/// Neighboring <see cref="ComputeNode"/>s of this <see cref="Node"/> with local vectors that have at least 1 common 
		/// entry with the local vector of this <see cref="Node"/>.
		/// </summary>
		internal SortedSet<int> ActiveNeighborsOfNode { get; }

		internal double[] InverseMultiplicities { get; }

		internal ComputeNode Node { get; }

		internal int NumIndices { get; }

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
			var activeNeighborsCopy = new SortedSet<int>(this.ActiveNeighborsOfNode);

			var inverseMultiplicitiesCopy = new double[this.InverseMultiplicities.Length];
			Array.Copy(this.InverseMultiplicities, inverseMultiplicitiesCopy, this.InverseMultiplicities.Length);

			var commonEntriesWithNeighborsCopy = new Dictionary<int, int[]>();
			foreach ((int nodeID, int[] data) in this.commonEntriesWithNeighbors)
			{
				var clonedData = new int[data.Length];
				Array.Copy(data, clonedData, data.Length);
				commonEntriesWithNeighborsCopy[nodeID] = clonedData;
			}

			return new LocalIndexer(
				this.Node,  activeNeighborsCopy, commonEntriesWithNeighborsCopy, inverseMultiplicitiesCopy, this.NumIndices);
		}


		internal int[] GetCommonEntriesWithNeighbor(int neighborID) => commonEntriesWithNeighbors[neighborID];

		private double[] FindMultiplicities()
		{
			var multiplicities = new int[NumIndices];
			for (int i = 0; i < NumIndices; ++i) multiplicities[i] = 1;
			foreach (int[] commonEntries in commonEntriesWithNeighbors.Values)
			{
				foreach (int i in commonEntries) multiplicities[i] += 1;
			}

			var inverseMultiplicities = new double[NumIndices];
			for (int i = 0; i < NumIndices; ++i) inverseMultiplicities[i] = 1.0 / multiplicities[i];

			return inverseMultiplicities;
		}
	}
}
