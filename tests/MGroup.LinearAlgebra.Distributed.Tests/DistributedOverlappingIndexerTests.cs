using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using MGroup.Environments;
using MGroup.LinearAlgebra.Distributed.Overlapping;
using MGroup.LinearAlgebra.Vectors;

using Xunit;

using static MGroup.LinearAlgebra.Distributed.Tests.Hexagon1DTestCase;

namespace MGroup.LinearAlgebra.Distributed.Tests
{
	public static class DistributedOverlappingIndexerTests
	{
		[Theory]
		[InlineData(EnvironmentChoice.SequentialSharedEnvironment)]
		[InlineData(EnvironmentChoice.TplSharedEnvironment)]
		public static void TestLocalNumEntriesManaged(EnvironmentChoice env) => TestLocalNumEntries(env.CreateEnvironment());

		internal static void TestLocalNumEntries(IComputeEnvironment environment)
		{
			environment.Initialize(CreateNodeTopology());
			DistributedOverlappingIndexer indexer = CreateIndexer(environment);

			environment.DoPerNode(nodeID =>
			{
				Assert.Equal(3, indexer.GetLocalComponent(nodeID).NumEntries);
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
				DistributedOverlappingIndexer.Local localIndexer = indexer.GetLocalComponent(nodeID);

				int[] neighborsExpected = GetNeighborsOfNode(nodeID);
				Assert.Equal(neighborsExpected, localIndexer.ActiveNeighborsOfNode.ToArray());
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
				DistributedOverlappingIndexer.Local localIndexer = indexer.GetLocalComponent(nodeID);

				Dictionary<int, int[]> commonEntriesExpected = CreateCommonEntriesWithNeighbors(nodeID);
				foreach (int neighborID in GetNeighborsOfNode(nodeID))
				{
					int[] commonEntriesComputed = localIndexer.GetCommonEntriesWithNeighbor(neighborID);
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
				DistributedOverlappingIndexer.Local localIndexer = indexer.GetLocalComponent(nodeID);

				var inverseMultiplicitiesExpected = Vector.CreateFromArray(new double[] { 0.5, 1, 0.5 });
				var inverseMultiplicitiesComputed = Vector.CreateFromArray(localIndexer.InverseMultiplicities);
				Assert.True(inverseMultiplicitiesExpected.Equals(inverseMultiplicitiesComputed, 1E-13));
			});
		}

		[InlineData(EnvironmentChoice.TplSharedEnvironment)]
		public static void TestCountCommonEntriesManaged(EnvironmentChoice env)
			=> TestCountCommonEntries(env.CreateEnvironment());

		internal static void TestCountCommonEntries(IComputeEnvironment environment)
		{
			environment.Initialize(CreateNodeTopology());
			DistributedOverlappingIndexer indexer = CreateIndexer(environment);

			environment.DoPerNode(nodeID =>
			{
				DistributedOverlappingIndexer.Local localIndexer = indexer.GetLocalComponent(nodeID);

				(int local, int remote) numCommonEntriesExpected = (1, 1);
				(int local, int remote) numCommonEntriesComputed = localIndexer.CountCommonEntries();
				Assert.Equal(numCommonEntriesExpected.local, numCommonEntriesComputed.local);
				Assert.Equal(numCommonEntriesExpected.remote, numCommonEntriesComputed.remote);
			});
		}

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
				DistributedOverlappingIndexer.Local localIndexer = indexer.GetLocalComponent(nodeID);

				ConcurrentDictionary<int, double[]> buffers = localIndexer.CreateBuffersForAllToAllWithNeighbors();
				int[] neighborIDs = GetNeighborsOfNode(nodeID);
				Assert.Equal(neighborIDs.Length, buffers.Count);

				foreach (int neighborID in neighborIDs)
				{
					int bufferSizeExpected = 1; // = num common entries between target node and ONE of its neighbors 
					Assert.Equal(bufferSizeExpected, buffers[neighborID].Length);
				}
			});
		}
	}
}
