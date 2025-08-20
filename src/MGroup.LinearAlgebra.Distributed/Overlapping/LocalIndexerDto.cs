using System;
using System.Collections.Generic;
using System.Text;

using MGroup.Environments;

namespace MGroup.LinearAlgebra.Distributed.Overlapping
{
	[Serializable]
    public class LocalIndexerDto
	{
		public LocalIndexerDto() { }

		internal LocalIndexerDto(LocalIndexer localIndexer)
		{
			NodeID = localIndexer.Node.ID;
			NumEntries = localIndexer.NumEntries;
			CommonEntriesOfNodeWithNeighbors = new Dictionary<int, int[]>();
			foreach (int neighborID in localIndexer.ActiveNeighborsOfNode)
			{
				this.CommonEntriesOfNodeWithNeighbors[neighborID] = localIndexer.GetCommonEntriesWithNeighbor(neighborID);
			}
		}

		public Dictionary<int, int[]> CommonEntriesOfNodeWithNeighbors { get; set; }

		public bool Modified { get; set; } = false;

		public int NodeID { get; set; }

		public int NumEntries { get; set; }

		internal LocalIndexer ToLocalIndexer(IComputeEnvironment environment)
		{
			var localIndexer = new LocalIndexer(environment.GetComputeNode(NodeID));
			localIndexer.Initialize(new LocalIndexerDto
			{
				NumEntries = this.NumEntries,
				CommonEntriesOfNodeWithNeighbors = this.CommonEntriesOfNodeWithNeighbors
			});
			return localIndexer;
		}
	}
}
