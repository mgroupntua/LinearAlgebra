using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;

using MGroup.Environments;
using MGroup.Environments.Mpi;
using MGroup.LinearAlgebra.Distributed.Overlapping;
using MGroup.LinearAlgebra.Iterative;
using MGroup.LinearAlgebra.Iterative.PreconditionedConjugateGradient;
using MGroup.LinearAlgebra.Iterative.PreconditionedConjugateGradient.BlockPcg;
using MGroup.LinearAlgebra.Iterative.Preconditioning;
using MGroup.LinearAlgebra.Iterative.Termination;
using MGroup.LinearAlgebra.Iterative.Termination.Iterations;
using MGroup.LinearAlgebra.Matrices;
using MGroup.LinearAlgebra.Vectors;

using Microsoft.VisualStudio.TestPlatform.PlatformAbstractions.Interfaces;

using Xunit;

using static MGroup.LinearAlgebra.Distributed.Tests.Overlapping.Hexagon1DTestCase;

namespace MGroup.LinearAlgebra.Distributed.Tests.Overlapping
{
	//TODO: Have multiple examples (IExample) that provide input and expected output for the methods that I want to test. Then 
	//      use theory with the examples and environments as params. Optionally the IExample will also provide a mock environment 
	//      different that has hardcoded logic that cannot fail and only works for the given example.
	public class DistributedOverlappingVectorTests
	{

		[Theory]
		[InlineData(EnvironmentChoice.SequentialSharedEnvironment)]
		[InlineData(EnvironmentChoice.TplSharedEnvironment)]
		public static void TestAxpyVectorsManaged(EnvironmentChoice env) => TestAxpyVectors(env.CreateEnvironment());

		internal static void TestAxpyVectors(IComputeEnvironment environment)
		{
			environment.Initialize(CreateNodeTopology());
			var indexer = CreateIndexer(environment);

			var localX = environment.CalcNodeData(n => GetX(n));
			var distributedX = new DistributedOverlappingVector(indexer, localX);

			var localY = environment.CalcNodeData(n => GetY(n));
			var distributedY = new DistributedOverlappingVector(indexer, localY);

			var localZExpected = environment.CalcNodeData(n => GetX(n) - 2.0 * GetY(n));
			var distributedZExpected = new DistributedOverlappingVector(indexer, localZExpected);

			var distributedZ = distributedX.CopyAsDistributed();
			distributedZ.AxpyIntoThis(distributedY, -2.0);

			var tol = 1E-13;
			Assert.True(distributedZExpected.Equals(distributedZ, tol));
		}

		[Theory]
		[InlineData(EnvironmentChoice.SequentialSharedEnvironment)]
		[InlineData(EnvironmentChoice.TplSharedEnvironment)]
		public static void TestDotProductManaged(EnvironmentChoice env) => TestDotProduct(env.CreateEnvironment());

		internal static void TestDotProduct(IComputeEnvironment environment)
		{
			environment.Initialize(CreateNodeTopology());
			var indexer = CreateIndexer(environment);

			var localX = environment.CalcNodeData(n => GetX(n));
			var distributedX = new DistributedOverlappingVector(indexer, localX);

			var localY = environment.CalcNodeData(n => GetY(n));
			var distributedY = new DistributedOverlappingVector(indexer, localY);

			var dotExpected = GetXDotY();
			var dot = distributedX.DotProduct(distributedY);

			var precision = 8;
			Assert.Equal(dotExpected, dot, precision);
		}

		[Theory]
		[InlineData(EnvironmentChoice.SequentialSharedEnvironment)]
		[InlineData(EnvironmentChoice.TplSharedEnvironment)]
		public static void TestEqualVectorsManaged(EnvironmentChoice env) => TestEqualVectors(env.CreateEnvironment());

		internal static void TestEqualVectors(IComputeEnvironment environment)
		{
			environment.Initialize(CreateNodeTopology());
			var indexer = CreateIndexer(environment);

			var localX = environment.CalcNodeData(n => GetX(n));
			var distributedX = new DistributedOverlappingVector(indexer, localX);

			var distributedXAlt = new DistributedOverlappingVector(indexer, localX);
			environment.DoPerNode(n => distributedXAlt.LocalVectors[n].CopyFrom(GetX(n)));

			var tol = 1E-13;
			Assert.True(distributedX.Equals(distributedXAlt, tol));
		}

		[Theory]
		[InlineData(EnvironmentChoice.SequentialSharedEnvironment)]
		[InlineData(EnvironmentChoice.TplSharedEnvironment)]
		public static void TestLinearCombinationVectorsManaged(EnvironmentChoice env)
			=> TestLinearCombinationVectors(env.CreateEnvironment());

		internal static void TestLinearCombinationVectors(IComputeEnvironment environment)
		{
			environment.Initialize(CreateNodeTopology());
			var indexer = CreateIndexer(environment);

			var localX = environment.CalcNodeData(n => GetX(n));
			var distributedX = new DistributedOverlappingVector(indexer, localX);

			var localY = environment.CalcNodeData(n => GetY(n));
			var distributedY = new DistributedOverlappingVector(indexer, localY);

			var localZExpected = environment.CalcNodeData(n => 2.0 * GetX(n) + 3.0 * GetY(n));
			var distributedZExpected = new DistributedOverlappingVector(indexer, localZExpected);

			var distributedZ = distributedX.CopyAsDistributed();
			distributedZ.LinearCombinationIntoThis(2.0, distributedY, 3.0);

			var tol = 1E-13;
			Assert.True(distributedZExpected.Equals(distributedZ, tol));
		}

		

		[Theory]
		[InlineData(EnvironmentChoice.SequentialSharedEnvironment)]
		[InlineData(EnvironmentChoice.TplSharedEnvironment)]
		public static void TestScaleVectorManaged(EnvironmentChoice env) => TestScaleVector(env.CreateEnvironment());

		internal static void TestScaleVector(IComputeEnvironment environment)
		{
			environment.Initialize(CreateNodeTopology());
			var indexer = CreateIndexer(environment);

			var localX = environment.CalcNodeData(n => GetX(n));
			var distributedX = new DistributedOverlappingVector(indexer, localX);

			var localZExpected = environment.CalcNodeData(n => -3.0 * GetX(n));
			var distributedZExpected = new DistributedOverlappingVector(indexer, localZExpected);

			var distributedZ = distributedX.CopyAsDistributed();
			distributedZ.ScaleIntoThis(-3.0);

			var tol = 1E-13;
			Assert.True(distributedZExpected.Equals(distributedZ, tol));
		}

		[Theory]
		[InlineData(EnvironmentChoice.SequentialSharedEnvironment)]
		[InlineData(EnvironmentChoice.TplSharedEnvironment)]
		public static void TestSumOverlappingEntriesManaged(EnvironmentChoice env) 
			=> TestSumOverlappingEntries(env.CreateEnvironment());

		internal static void TestSumOverlappingEntries(IComputeEnvironment environment)
		{
			environment.Initialize(CreateNodeTopology());
			var indexer = CreateIndexer(environment);

			var localInputW = environment.CalcNodeData(n => GetWBeforeSumOverlapping(n));
			var distributedInputW = new DistributedOverlappingVector(indexer, localInputW);

			var localOutputW = environment.CalcNodeData(n => GetWAfterSumOverlapping(n));
			var distributedOutputW = new DistributedOverlappingVector(indexer, localOutputW);

			distributedInputW.SumOverlappingEntries();

			var tol = 1E-13;
			Assert.True(distributedOutputW.Equals(distributedInputW, tol));
		}

		public static void RunMpiTests()
		{
			// Launch 3 processes
			using (var mpiEnvironment = new MpiEnvironment(new MasterSlavesGlobalOperationStrategy()))
			{
				MpiDebugUtilities.AssistDebuggerAttachment();

				TestAxpyVectors(mpiEnvironment);
				TestDotProduct(mpiEnvironment);
				TestEqualVectors(mpiEnvironment);
				TestLinearCombinationVectors(mpiEnvironment);
				TestScaleVector(mpiEnvironment);
				TestSumOverlappingEntries(mpiEnvironment);

				MpiDebugUtilities.DoSerially(MPI.Communicator.world,
					() => Console.WriteLine($"Process {MPI.Communicator.world.Rank}: All tests passed"));
			}
		}

		internal static Vector GetX(int nodeID)
		{
			var x = Vector.CreateZero(3);
			for (var i = 0; i < x.Length; ++i)
			{
				x[i] = 2 * nodeID + i;
			}
			if (nodeID == NumComputeNodes - 1)
			{
				x[x.Length - 1] = 0;
			}
			return x;
		}

		private static Vector GetY(int nodeID)
		{
			var y = GetX(nodeID);
			y.DoToAllEntriesIntoThis(a => 10 + a);
			return y;
		}

		private static double GetXDotY() => 1166;

		private static Vector GetWBeforeSumOverlapping(int nodeID)
		{
			var w = Vector.CreateZero(3);
			for (var i = 0; i < w.Length; ++i)
			{
				w[i] = 3 * nodeID + i;
			}
			return w;
		}

		private static Vector GetWAfterSumOverlapping(int nodeID)
		{
			var previous = nodeID - 1 >= 0 ? nodeID - 1 : NumComputeNodes - 1;
			var next = (nodeID + 1) % NumComputeNodes;

			var w = GetWBeforeSumOverlapping(nodeID);
			var wPrevious = GetWBeforeSumOverlapping(previous);
			var wNext = GetWBeforeSumOverlapping(next);

			w[0] += wPrevious[2];
			w[2] += wNext[0];

			return w;
		}
	}
}
