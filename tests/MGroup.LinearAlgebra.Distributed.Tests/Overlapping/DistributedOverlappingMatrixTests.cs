using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using MGroup.Environments;
using MGroup.Environments.Mpi;
using MGroup.LinearAlgebra.Distributed.Overlapping;
using MGroup.LinearAlgebra.Matrices;
using MGroup.LinearAlgebra.Vectors;

using Xunit;

using static MGroup.LinearAlgebra.Distributed.Tests.Overlapping.Hexagon1DTestCase;

namespace MGroup.LinearAlgebra.Distributed.Tests.Overlapping
{
    public class DistributedOverlappingMatrixTests
    {
		[Theory]
		[InlineData(EnvironmentChoice.SequentialSharedEnvironment)]
		[InlineData(EnvironmentChoice.TplSharedEnvironment)]
		public static void TestMatrixVectorMultiplicationManaged(EnvironmentChoice env)
			=> TestMatrixVectorMultiplication(env.CreateEnvironment());

		internal static void TestMatrixVectorMultiplication(IComputeEnvironment environment)
		{
			environment.Initialize(CreateNodeTopology());
			DistributedOverlappingIndexer indexer = CreateIndexer(environment);

			var distributedA = new DistributedOverlappingTransformation(indexer,
				(n, x, y) => GetMatrixA(n).MultiplyIntoResult(x, y));

			Dictionary<int, Vector> localX = environment.CalcNodeData(n => GetVectorX(n));
			var distributedX = new DistributedOverlappingVector(indexer, localX);

			Dictionary<int, Vector> localAxExpected = environment.CalcNodeData(n => GetVectorAx(n));
			var distributedAxExpected = new DistributedOverlappingVector(indexer, localAxExpected);

			var distributedAx = new DistributedOverlappingVector(indexer);
			distributedA.Multiply(distributedX, distributedAx);

			var tol = 1E-13;
			Assert.True(distributedAxExpected.Equals(distributedAx, tol));
		}

		public static void RunMpiTests()
		{
			// Launch 3 processes
			using (var mpiEnvironment = new MpiEnvironment(new MasterSlavesGlobalOperationStrategy()))
			{
				MpiDebugUtilities.AssistDebuggerAttachment();

				TestMatrixVectorMultiplication(mpiEnvironment);

				MpiDebugUtilities.DoSerially(MPI.Communicator.world,
					() => Console.WriteLine($"Process {MPI.Communicator.world.Rank}: All tests passed"));
			}
		}

		internal static Matrix GetMatrixA(int nodeID)
		{
			double[,] A;
			if (nodeID == 0)
			{
				A = new double[,]
				{
					{ 100, 1, 2 },
					{ 1, 101, 3 },
					{ 2, 3, 102 }
				};
			}
			else if (nodeID == 1)
			{
				A = new double[,]
				{
					{ 103, 6, 7 },
					{ 6, 104, 8 },
					{ 7, 8, 105 }
				};
			}
			else if (nodeID == 2)
			{
				A = new double[,]
				{
					{ 106, 11, 12 },
					{ 11, 107, 13 },
					{ 12, 13, 108 }
				};
			}
			else if (nodeID == 3)
			{
				A = new double[,]
				{
					{ 109, 16, 17 },
					{ 16, 110, 18 },
					{ 17, 18, 111 }
				};
			}
			else if (nodeID == 4)
			{
				A = new double[,]
				{
					{ 112, 21, 22 },
					{ 21, 113, 23 },
					{ 22, 23, 114 }
				};
			}
			else
			{
				Debug.Assert(nodeID == 5);
				A = new double[,]
				{
					{ 115, 26, 15 },
					{ 26, 116, 16 },
					{ 15, 16, 105 }
				};
			}
			return Matrix.CreateFromArray(A);
		}

		internal static Vector GetVectorAx(int nodeID)
		{
			double[] Ax;
			if (nodeID == 0)
			{
				Ax = new double[] { 331, 107, 459 };
			}
			else if (nodeID == 1)
			{
				Ax = new double[] { 459, 356, 1009 };
			}
			else if (nodeID == 2)
			{
				Ax = new double[] { 1009, 657, 1663 };
			}
			else if (nodeID == 3)
			{
				Ax = new double[] { 1663, 1010, 2421 };
			}
			else if (nodeID == 4)
			{
				Ax = new double[] { 2421, 1415, 2959 };
			}
			else
			{
				Debug.Assert(nodeID == 5);
				Ax = new double[] { 2959, 1536, 331 };
			}
			return Vector.CreateFromArray(Ax);
		}

		private static Vector GetVectorX(int nodeID) => DistributedOverlappingVectorTests.GetX(nodeID);
	}
}
