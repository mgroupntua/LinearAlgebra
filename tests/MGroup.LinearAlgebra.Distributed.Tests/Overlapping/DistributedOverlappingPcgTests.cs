using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using MGroup.Environments;
using MGroup.LinearAlgebra.Distributed.Overlapping;
using MGroup.LinearAlgebra.Iterative.PreconditionedConjugateGradient.BlockPcg;
using MGroup.LinearAlgebra.Iterative.PreconditionedConjugateGradient;
using MGroup.LinearAlgebra.Iterative.Preconditioning;
using MGroup.LinearAlgebra.Iterative.Termination.Iterations;
using MGroup.LinearAlgebra.Iterative;
using Xunit;
using MGroup.Environments.Mpi;

using static MGroup.LinearAlgebra.Distributed.Tests.Overlapping.Hexagon1DTestCase;
using MGroup.LinearAlgebra.Vectors;
using MGroup.LinearAlgebra.Matrices;

namespace MGroup.LinearAlgebra.Distributed.Tests.Overlapping
{
    public class DistributedOverlappingPcgTests
    {
		[Theory]
		[InlineData(EnvironmentChoice.SequentialSharedEnvironment)]
		[InlineData(EnvironmentChoice.TplSharedEnvironment)]
		public static void TestPcgManaged(EnvironmentChoice env) => TestPcg(env.CreateEnvironment());

		internal static void TestPcg(IComputeEnvironment environment)
		{
			environment.Initialize(CreateNodeTopology());
			DistributedOverlappingIndexer indexer = CreateIndexer(environment);

			var distributedA = new DistributedOverlappingTransformation(indexer,
				(n, x, y) => GetMatrixA(n).MultiplyIntoResult(x, y));

			Dictionary<int, Vector> localAx = environment.CalcNodeData(n => GetVectorAx(n));
			var distributedAx = new DistributedOverlappingVector(indexer, localAx);

			Dictionary<int, Vector> localXExpected = environment.CalcNodeData(n => GetVectorX(n));
			var distributedXExpected = new DistributedOverlappingVector(indexer, localXExpected);

			var pcgFactory = new PcgAlgorithm.Factory();
			var maxIterations = 12;
			pcgFactory.MaxIterationsProvider = new FixedMaxIterationsProvider(maxIterations);
			pcgFactory.ResidualTolerance = 1E-10;
			var pcg = pcgFactory.Build();
			var distributedX = new DistributedOverlappingVector(indexer);
			var stats = pcg.Solve(distributedA, new IdentityPreconditioner(), distributedAx, distributedX, true);

			var tol = 1E-10;
			Assert.True(stats.HasConverged);
			Assert.True(distributedXExpected.Equals(distributedX, tol));
		}

		[Theory]
		[InlineData(EnvironmentChoice.SequentialSharedEnvironment)]
		[InlineData(EnvironmentChoice.TplSharedEnvironment)]
		public static void TestBlockPcgManaged(EnvironmentChoice env) => TestBlockPcg(env.CreateEnvironment());

		internal static void TestBlockPcg(IComputeEnvironment environment)
		{
			environment.Initialize(CreateNodeTopology());
			DistributedOverlappingIndexer indexer = CreateIndexer(environment);

			var distributedA = new DistributedOverlappingTransformation(indexer,
				(n, x, y) => GetMatrixA(n).MultiplyIntoResult(x, y));

			Dictionary<int, Vector> localAx = environment.CalcNodeData(n => GetVectorAx(n));
			var distributedAx = new DistributedOverlappingVector(indexer, localAx);

			Dictionary<int, Vector> localXExpected = environment.CalcNodeData(n => GetVectorX(n));
			var distributedXExpected = new DistributedOverlappingVector(indexer, localXExpected);

			var pcgFactory = new BlockPcgAlgorithm.Factory();
			var maxIterations = 12;
			pcgFactory.MaxIterationsProvider = new FixedMaxIterationsProvider(maxIterations);
			pcgFactory.ResidualTolerance = 1E-10;
			var pcg = pcgFactory.Build();
			var distributedX = new DistributedOverlappingVector(indexer);
			var stats = pcg.Solve(distributedA, new IdentityPreconditioner(), distributedAx, distributedX, true);

			var tol = 1E-6;
			Assert.True(stats.HasConverged);
			Assert.True(distributedXExpected.Equals(distributedX, tol));
		}

		public static void RunMpiTests()
		{
			// Launch 3 processes
			using (var mpiEnvironment = new MpiEnvironment(new MasterSlavesGlobalOperationStrategy()))
			{
				MpiDebugUtilities.AssistDebuggerAttachment();

				TestBlockPcg(mpiEnvironment);
				TestPcg(mpiEnvironment);

				MpiDebugUtilities.DoSerially(MPI.Communicator.world,
					() => Console.WriteLine($"Process {MPI.Communicator.world.Rank}: All tests passed"));
			}
		}

		private static Matrix GetMatrixA(int nodeID) => DistributedOverlappingMatrixTests.GetMatrixA(nodeID);

		private static Vector GetVectorAx(int nodeID) => DistributedOverlappingMatrixTests.GetVectorAx(nodeID);

		private static Vector GetVectorX(int nodeID) => DistributedOverlappingVectorTests.GetX(nodeID);



	}
}
