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
// Dofs: 1 per node
namespace MGroup.LinearAlgebra.Distributed.Tests
{
    public static class Hexagon1DTestCase
    {
		public static int NumNodes => 6;

		public static int[,] AllNodeNeighbors => new int[,]
		{
			{ 5, 1},
			{ 0, 2 },
			{ 1, 3 },
			{ 2, 4 },
			{ 3, 5 },
			{ 4, 0 }
		};

		public static Dictionary<int, int[]> CreateCommonEntriesWithNeighbors(int nodeID)
		{
			int previous = nodeID - 1 >= 0 ? nodeID - 1 : NumNodes - 1;
			int next = (nodeID + 1) % NumNodes;
			var commonEntries = new Dictionary<int, int[]>();
			commonEntries[previous] = new int[] { 0 };
			commonEntries[next] = new int[] { 2 };
			return commonEntries;
		}

		public static DistributedOverlappingIndexer CreateIndexer(IComputeEnvironment environment)
		{
			var indexer = new DistributedOverlappingIndexer(environment);
			Action<int> initializeIndexer = n =>
			{
				Dictionary<int, int[]> commonEntries = CreateCommonEntriesWithNeighbors(n);
				int numEntries = 3;
				indexer.GetLocalComponent(n).Initialize(numEntries, commonEntries);
			};
			environment.DoPerNode(initializeIndexer);

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
