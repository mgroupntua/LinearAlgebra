using System;
using System.Collections.Generic;
using System.Text;

using MGroup.Environments;

namespace MGroup.LinearAlgebra.Distributed.Overlapping
{
	[Serializable]
    public class LocalIndexerDto
	{
		internal LocalIndexerDto() { }
		public static LocalIndexerDto CreateUnmodified() => new LocalIndexerDto { Modified = false };

		public static LocalIndexerDto CreateWithNewContent(
			int numIndices, Dictionary<int, int[]> commonEntriesOfNodeWithNeighbors)
		{
			var dto = new LocalIndexerDto
			{
				NumIndices = numIndices,
				CommonEntriesOfNodeWithNeighbors = commonEntriesOfNodeWithNeighbors
			};

			return dto;
		}

		internal static LocalIndexerDto CreateForSerialization(LocalIndexer localIndexer)
		{
			var dto = new LocalIndexerDto
			{
				NodeID = localIndexer.Node.ID,
				NumIndices = localIndexer.NumIndices,
				CommonEntriesOfNodeWithNeighbors = new Dictionary<int, int[]>()
			};

			foreach (int neighborID in localIndexer.ActiveNeighborsOfNode)
			{
				dto.CommonEntriesOfNodeWithNeighbors[neighborID] = localIndexer.GetCommonEntriesWithNeighbor(neighborID);
			}

			return dto;
		}

		public Dictionary<int, int[]> CommonEntriesOfNodeWithNeighbors { get; set; }

		public bool Modified { get; set; } = true;

		public int NodeID { get; set; } = -1;

		public int NumIndices { get; set; }

		internal LocalIndexer ToLocalIndexer(IComputeEnvironment environment)
		{
			ComputeNode node = environment.GetComputeNode(this.NodeID);
			return new LocalIndexer(node, this.CommonEntriesOfNodeWithNeighbors, this.NumIndices);
		}
	}
}
