using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using MGroup.Environments;
using MGroup.LinearAlgebra.Distributed.Overlapping;

//          8          
//        /   \          
//      9       7       
//     /         \       
//   10           6                          
//    |           |      
//   11           5     s0:        s1:        s2:        s3:        s4:        s5:      
//    |           |      
//    0           4     0                4     6         8                 8   10
//     \         /       \              /      |           \             /      |
//       1     3           1          3        5             7         9       11
//        \   /             \        /         |              \       /         |
//          2                 2    2           4               6    10          0
//
// Clusters:  c0 (s0, s1), c1 (s2, s3), c2 (s4, s5) 
// Dofs: 1 per physical node
namespace MGroup.LinearAlgebra.Distributed.Tests.Overlapping
{
    public static class Hexagon1DTestCase
    {
		public static int NumComputeNodes => 6;

		public static int[,] AllNodeNeighbors => new int[,]
		{
			{ 5, 1 },
			{ 0, 2 },
			{ 1, 3 },
			{ 2, 4 },
			{ 3, 5 },
			{ 4, 0 }
		};

		public static List<Dictionary<int, int>> GlobalToLocalIndices 
		{
			get
			{
				var result = new List<Dictionary<int, int>>();

				for (var g = 0; g < 12; g++)
				{
					result.Add(new Dictionary<int, int>());
				}

				// s0
				result[0].Add(0, 0);
				result[1].Add(0, 1);
				result[2].Add(0, 2);

				// s1
				result[2].Add(1, 0);
				result[3].Add(1, 1);
				result[4].Add(1, 2);

				// s2
				result[4].Add(2, 0);
				result[5].Add(2, 1);
				result[6].Add(2, 2);

				// s3
				result[6].Add(3, 0);
				result[7].Add(3, 1);
				result[8].Add(3, 2);

				// s4
				result[8].Add(4, 0);
				result[9].Add(4, 1);
				result[10].Add(4, 2);

				// s5
				result[10].Add(5, 0);
				result[11].Add(5, 1);
				result[0].Add(5, 2);

				return result;
			}
		}

		/// <summary>
		/// Row = compute node ID. Column = local index. Result = global index
		/// </summary>
		public static int[][] LocalToGlobalIndices => new int[][]
		{
			new int[] { 0, 1, 2 },
			new int[] { 2, 3, 4 },
			new int[] { 4, 5, 6 },
			new int[] { 6, 7, 8 },
			new int[] { 8, 9, 10 },
			new int[] { 10, 11, 0 }
		};

		public static Dictionary<int, int[]> CreateCommonEntriesWithNeighbors(int nodeID)
		{
			var previous = nodeID - 1 >= 0 ? nodeID - 1 : NumComputeNodes - 1;
			var next = (nodeID + 1) % NumComputeNodes;
			var commonEntries = new Dictionary<int, int[]>();
			commonEntries[previous] = new int[] { 0 };
			commonEntries[next] = new int[] { 2 };
			return commonEntries;
		}

		public static DistributedOverlappingIndexer CreateIndexer(IComputeEnvironment environment)
		{
			var indexer = new DistributedOverlappingIndexer(environment);
			indexer.Initialize(nodeID => new LocalIndexerDto
			{
				NumEntries = 3,
				CommonEntriesOfNodeWithNeighbors = CreateCommonEntriesWithNeighbors(nodeID),
			});
			return indexer;
		}


		public static int[] GetNeighborsOfNode(int nodeID)
		{
			var result = new int[AllNodeNeighbors.GetLength(1)];
			for (var i = 0; i < result.Length; i++)
			{
				result[i] = AllNodeNeighbors[nodeID, i];
			}

			Array.Sort(result);
			return result;
		}

		public static ComputeNodeTopology CreateNodeTopology() //TODOMPI: this is also described in Environments.Tests. 
		{
			var topology = new ComputeNodeTopology();
			topology.AddNode(0, new int[] { 5, 1 }, 0);
			topology.AddNode(1, new int[] { 0, 2 }, 0);
			topology.AddNode(2, new int[] { 1, 3 }, 1);
			topology.AddNode(3, new int[] { 2, 4 }, 1);
			topology.AddNode(4, new int[] { 3, 5 }, 2);
			topology.AddNode(5, new int[] { 4, 0 }, 2);
			return topology;
		}
	}
}
