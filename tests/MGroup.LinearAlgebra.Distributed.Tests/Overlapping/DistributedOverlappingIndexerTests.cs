using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using MGroup.Environments;
using MGroup.Environments.Mpi;
using MGroup.LinearAlgebra.Distributed.Overlapping;
using MGroup.LinearAlgebra.Vectors;

using Xunit;

using static MGroup.LinearAlgebra.Distributed.Overlapping.DistributedOverlappingIndexer;
using static MGroup.LinearAlgebra.Distributed.Tests.Overlapping.Hexagon1DTestCase;

namespace MGroup.LinearAlgebra.Distributed.Tests.Overlapping
{
	public static class DistributedOverlappingIndexerTests
	{
		[Theory]
		[InlineData(EnvironmentChoice.SequentialSharedEnvironment)]
		[InlineData(EnvironmentChoice.TplSharedEnvironment)]
		public static void TestCreateAllToAllBuffersManaged(EnvironmentChoice env)
			=> TestCreateAllToAllBuffers(env.CreateEnvironment());

		internal static void TestCreateAllToAllBuffers(IComputeEnvironment environment)
		{
			environment.Initialize(CreateNodeTopology());
			DistributedOverlappingIndexer indexer = CreateIndexer(environment);

			environment.DoPerNode(nodeID =>
			{
				var buffers = indexer.CreateBuffersForAllToAllWithNeighbors(nodeID);
				var neighborIDs = GetNeighborsOfNode(nodeID);
				Assert.Equal(neighborIDs.Length, buffers.Count);

				foreach (var neighborID in neighborIDs)
				{
					var bufferSizeExpected = 1; // = num common entries between target node and ONE of its neighbors 
					Assert.Equal(bufferSizeExpected, buffers[neighborID].Length);
				}
			});
		}

		[InlineData(EnvironmentChoice.SequentialSharedEnvironment)]
		[InlineData(EnvironmentChoice.TplSharedEnvironment)]
		public static void TestCountCommonEntriesManaged(EnvironmentChoice env)
			=> TestCountCommonEntries(env.CreateEnvironment());

		internal static void TestCountCommonEntries(IComputeEnvironment environment)
		{
			environment.Initialize(CreateNodeTopology());
			DistributedOverlappingIndexer indexer = CreateIndexer(environment);

			environment.DoPerNode(nodeID =>
			{
				(int local, int remote) numCommonEntriesExpected = (1, 1);
				var numCommonEntriesComputed = indexer.CountCommonEntriesOfNodeWithNeighbors(nodeID);
				Assert.Equal(numCommonEntriesExpected.local, numCommonEntriesComputed.local);
				Assert.Equal(numCommonEntriesExpected.remote, numCommonEntriesComputed.remote);
			});
		}

		[Theory]
		[InlineData(EnvironmentChoice.SequentialSharedEnvironment)]
		[InlineData(EnvironmentChoice.TplSharedEnvironment)]
		public static void TestGlobalToLocalIndexManaged(EnvironmentChoice env)
			=> TestGlobalToLocalIndex(env.CreateEnvironment());

		internal static void TestGlobalToLocalIndex(IComputeEnvironment environment)
		{
			environment.Initialize(CreateNodeTopology());
			var indexer = CreateIndexer(environment);

			for (var globalIdx = 0; globalIdx < indexer.NumGlobalIndices; globalIdx++)
			{
				for (var node = 0; node <= NumComputeNodes; node++)
				{
					var localIdxComputed = indexer.FindLocalIndexOf(globalIdx, node);

					var localIdxExpected = -1;
					if (GlobalToLocalIndices[globalIdx].ContainsKey(node))
					{
						localIdxExpected = GlobalToLocalIndices[globalIdx][node];
					}

					Assert.Equal(localIdxExpected, localIdxComputed);
				}
			}
		}

		[Theory]
		[InlineData(EnvironmentChoice.SequentialSharedEnvironment)]
		[InlineData(EnvironmentChoice.TplSharedEnvironment)]
		public static void TestGlobalToLocalIndicesManaged(EnvironmentChoice env)
			=> TestGlobalToLocalIndices(env.CreateEnvironment());

		internal static void TestGlobalToLocalIndices(IComputeEnvironment environment)
		{
			environment.Initialize(CreateNodeTopology());
			var indexer = CreateIndexer(environment);

			for (var gi = 0; gi < indexer.NumGlobalIndices; gi++)
			{
				var localIndicesComputed = indexer.FindLocalIndicesOf(gi);
				var localIndicesExpected = GlobalToLocalIndices[gi];
				Assert.Equal(localIndicesExpected.Count, localIndicesComputed.Count);
				foreach ((var nodeID, var localIdx) in localIndicesComputed)
				{
					Assert.True(localIndicesExpected.ContainsKey(nodeID));
					Assert.Equal(localIndicesExpected[nodeID], localIdx);
				}
			}
		}

		[Theory]
		[InlineData(EnvironmentChoice.SequentialSharedEnvironment)]
		[InlineData(EnvironmentChoice.TplSharedEnvironment)]
		public static void TestLocalToGlobalIndexManaged(EnvironmentChoice env) 
			=> TestLocalToGlobalIndex(env.CreateEnvironment());

		internal static void TestLocalToGlobalIndex(IComputeEnvironment environment)
		{
			environment.Initialize(CreateNodeTopology());
			var indexer = CreateIndexer(environment);

			environment.DoPerNode(nodeID =>
			{
				var numLocalIndices = indexer.GetNumLocalIndices(nodeID);
				var localIndices = Enumerable.Range(0, numLocalIndices).ToArray();
				var globalIndicesComputed = localIndices.Select(li => indexer.FindGlobalIndexOf(nodeID, li)).ToArray();
				var globalIndicesExpected = localIndices.Select(li => LocalToGlobalIndices[nodeID][li]).ToArray();
				Assert.Equal(globalIndicesExpected, globalIndicesComputed);
			});
		}

		[Theory]
		[InlineData(EnvironmentChoice.SequentialSharedEnvironment)]
		[InlineData(EnvironmentChoice.TplSharedEnvironment)]
		public static void TestLocalNeighborsManaged(EnvironmentChoice env) => TestLocalNeighbors(env.CreateEnvironment());

		internal static void TestLocalNeighbors(IComputeEnvironment environment)
		{
			environment.Initialize(CreateNodeTopology());
			DistributedOverlappingIndexer indexer = CreateIndexer(environment);

			environment.DoPerNode(nodeID =>
			{
				var neighborsExpected = GetNeighborsOfNode(nodeID);
				Assert.Equal(neighborsExpected, indexer.GetActiveNeighborIDs(nodeID).ToArray());
			});
		}

		[Theory]
		[InlineData(EnvironmentChoice.SequentialSharedEnvironment)]
		[InlineData(EnvironmentChoice.TplSharedEnvironment)]
		public static void TestLocalCommonEntriesManaged(EnvironmentChoice env) 
			=> TestLocalCommonEntries(env.CreateEnvironment());

		internal static void TestLocalCommonEntries(IComputeEnvironment environment)
		{
			environment.Initialize(CreateNodeTopology());
			DistributedOverlappingIndexer indexer = CreateIndexer(environment);

			environment.DoPerNode(nodeID =>
			{
				Dictionary<int, int[]> commonEntriesExpected = CreateCommonEntriesWithNeighbors(nodeID);
				foreach (var neighborID in GetNeighborsOfNode(nodeID))
				{
					var commonEntriesComputed = indexer.GetCommonEntriesOfNodeWithNeighbor(nodeID, neighborID);
					Assert.Equal(commonEntriesExpected[neighborID], commonEntriesComputed);
				}
			});
		}

		[Theory]
		[InlineData(EnvironmentChoice.SequentialSharedEnvironment)]
		[InlineData(EnvironmentChoice.TplSharedEnvironment)]
		public static void TestLocalMultiplicitiesManaged(EnvironmentChoice env)
			=> TestLocalMultiplicities(env.CreateEnvironment());

		internal static void TestLocalMultiplicities(IComputeEnvironment environment)
		{
			environment.Initialize(CreateNodeTopology());
			DistributedOverlappingIndexer indexer = CreateIndexer(environment);

			environment.DoPerNode(nodeID =>
			{
				var inverseMultiplicitiesExpected = Vector.CreateFromArray(new double[] { 0.5, 1, 0.5 });
				var inverseMultiplicitiesComputed = Vector.CreateFromArray(indexer.GetInverseMultiplicities(nodeID));
				Assert.True(inverseMultiplicitiesExpected.Equals(inverseMultiplicitiesComputed, 1E-13));
			});
		}

		[Theory]
		[InlineData(EnvironmentChoice.SequentialSharedEnvironment)]
		[InlineData(EnvironmentChoice.TplSharedEnvironment)]
		public static void TestNumGlobalEntriesManaged(EnvironmentChoice env) => TestNumGlobalEntries(env.CreateEnvironment());

		internal static void TestNumGlobalEntries(IComputeEnvironment environment)
		{
			environment.Initialize(CreateNodeTopology());
			DistributedOverlappingIndexer indexer = CreateIndexer(environment);
			Assert.Equal(NumGlobalEntries, indexer.NumGlobalIndices);
			
		}

		[Theory]
		[InlineData(EnvironmentChoice.SequentialSharedEnvironment)]
		[InlineData(EnvironmentChoice.TplSharedEnvironment)]
		public static void TestNumLocalEntriesManaged(EnvironmentChoice env) => TestNumLocalEntries(env.CreateEnvironment());

		internal static void TestNumLocalEntries(IComputeEnvironment environment)
		{
			environment.Initialize(CreateNodeTopology());
			DistributedOverlappingIndexer indexer = CreateIndexer(environment);

			environment.DoPerNode(nodeID =>
			{
				Assert.Equal(3, indexer.GetNumLocalIndices(nodeID));
			});
		}

		[Theory]
		[InlineData(EnvironmentChoice.SequentialSharedEnvironment)]
		[InlineData(EnvironmentChoice.TplSharedEnvironment)]
		public static void TestReuseAsBasisForNewIndexerManaged(EnvironmentChoice env)
			=> TestReuseAsBasisForNewIndexer(env.CreateEnvironment());

		internal static void TestReuseAsBasisForNewIndexer(IComputeEnvironment environment)
		{
			environment.Initialize(CreateNodeTopology());
			DistributedOverlappingIndexer originalIndexer = CreateIndexer(environment);

			// Add an internal entry to the last subdomain
			DistributedOverlappingIndexer newIndexer = originalIndexer.ReuseAsBasisForNewIndexer(nodeID =>
			{
				if (nodeID == 5) 
				{
					Dictionary<int, int[]> commonEntriesWithNodes = [];
					commonEntriesWithNodes[4] = new int[] { 0 };
					commonEntriesWithNodes[0] = new int[] { 3 };
					return LocalIndexerDto.CreateWithNewContent(4, commonEntriesWithNodes);
				}
				else
				{
					return LocalIndexerDto.CreateUnmodified();
				}
			});

			environment.DoPerNode(nodeID =>
			{
				var inverseMultiplicitiesExpected = Vector.CreateFromArray(
					nodeID == 5 ? new double[] { 0.5, 1, 1, 0.5 } : new double[] { 0.5, 1, 0.5 });
				var inverseMultiplicitiesComputed = Vector.CreateFromArray(newIndexer.GetInverseMultiplicities(nodeID));
				Assert.True(inverseMultiplicitiesExpected.Equals(inverseMultiplicitiesComputed, 1E-13));
			});
			Assert.Equal(NumGlobalEntries + 1, newIndexer.NumGlobalIndices);
		}

		public static void RunMpiTests()
		{
			// Launch 3 processes
			using (var mpiEnvironment = new MpiEnvironment(new MasterSlavesGlobalOperationStrategy()))
			{
				MpiDebugUtilities.AssistDebuggerAttachment();

				TestCountCommonEntries(mpiEnvironment);
				TestCreateAllToAllBuffers(mpiEnvironment);
				TestLocalCommonEntries(mpiEnvironment);
				TestLocalNeighbors(mpiEnvironment);
				TestNumLocalEntries(mpiEnvironment);
				TestLocalMultiplicities(mpiEnvironment);
				TestGlobalToLocalIndices(mpiEnvironment);
				TestLocalToGlobalIndex(mpiEnvironment);
				TestReuseAsBasisForNewIndexer(mpiEnvironment);

				MpiDebugUtilities.DoSerially(MPI.Communicator.world,
					() => Console.WriteLine($"Process {MPI.Communicator.world.Rank}: All tests passed"));
			}
		}
	}
}
