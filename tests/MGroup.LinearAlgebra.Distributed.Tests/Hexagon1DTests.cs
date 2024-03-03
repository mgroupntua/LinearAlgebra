using System;
using System.Collections.Generic;
using System.Diagnostics;

using MGroup.Environments;
using MGroup.LinearAlgebra.Distributed.IterativeMethods.PCG;
using MGroup.LinearAlgebra.Distributed.IterativeMethods.Preconditioning;
using MGroup.LinearAlgebra.Distributed.Overlapping;
using MGroup.LinearAlgebra.Iterative;
using MGroup.LinearAlgebra.Iterative.Termination;
using MGroup.LinearAlgebra.Iterative.Termination.Iterations;
using MGroup.LinearAlgebra.Matrices;
using MGroup.LinearAlgebra.Vectors;

using Xunit;

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
// 1 dof per node
namespace MGroup.LinearAlgebra.Distributed.Tests
{
	//TODO: Have multiple examples (IExample) that provide input and expected output for the methods that I want to test. Then 
	//      use theory with the examples and environments as params. Optionally the IExample will also provide a mock environment 
	//      different that has hardcoded logic that cannot fail and only works for the given example.
	public class Hexagon1DTests
	{
		public const int numNodes = 6;

		[Theory]
		[InlineData(EnvironmentChoice.SequentialSharedEnvironment)]
		[InlineData(EnvironmentChoice.TplSharedEnvironment)]
		public static void TestAxpyVectorsManaged(EnvironmentChoice env) => TestAxpyVectors(env.CreateEnvironment());

		internal static void TestAxpyVectors(IComputeEnvironment environment)
		{
			environment.Initialize(CreateNodeTopology());
			DistributedOverlappingIndexer indexer = CreateIndexer(environment);

			Dictionary<int, Vector> localX = environment.CalcNodeData(n => GetX(n));
			var distributedX = new DistributedOverlappingVector(indexer, localX);

			Dictionary<int, Vector> localY = environment.CalcNodeData(n => GetY(n));
			var distributedY = new DistributedOverlappingVector(indexer, localY);

			Dictionary<int, Vector> localZExpected = environment.CalcNodeData(n => GetX(n) - 2.0 * GetY(n));
			var distributedZExpected = new DistributedOverlappingVector(indexer, localZExpected);

			DistributedOverlappingVector distributedZ = distributedX.Copy();
			distributedZ.AxpyIntoThis(distributedY, -2.0);

			double tol = 1E-13;
			Assert.True(distributedZExpected.Equals(distributedZ, tol));
		}

		[Theory]
		[InlineData(EnvironmentChoice.SequentialSharedEnvironment)]
		[InlineData(EnvironmentChoice.TplSharedEnvironment)]
		public static void TestDotProductManaged(EnvironmentChoice env) => TestDotProduct(env.CreateEnvironment());

		internal static void TestDotProduct(IComputeEnvironment environment)
		{
			environment.Initialize(CreateNodeTopology());
			DistributedOverlappingIndexer indexer = CreateIndexer(environment);

			Dictionary<int, Vector> localX = environment.CalcNodeData(n => GetX(n));
			var distributedX = new DistributedOverlappingVector(indexer, localX);

			Dictionary<int, Vector> localY = environment.CalcNodeData(n => GetY(n));
			var distributedY = new DistributedOverlappingVector(indexer, localY);

			double dotExpected = GetXDotY();
			double dot = distributedX.DotProduct(distributedY);

			int precision = 8;
			Assert.Equal(dotExpected, dot, precision);
		}

		[Theory]
		[InlineData(EnvironmentChoice.SequentialSharedEnvironment)]
		[InlineData(EnvironmentChoice.TplSharedEnvironment)]
		public static void TestEqualVectorsManaged(EnvironmentChoice env) => TestEqualVectors(env.CreateEnvironment());

		internal static void TestEqualVectors(IComputeEnvironment environment)
		{
			environment.Initialize(CreateNodeTopology());
			DistributedOverlappingIndexer indexer = CreateIndexer(environment);

			Dictionary<int, Vector> localX = environment.CalcNodeData(n => GetX(n));
			var distributedX = new DistributedOverlappingVector(indexer, localX);

			var distributedXAlt = new DistributedOverlappingVector(indexer, localX);
			environment.DoPerNode(n => distributedXAlt.LocalVectors[n].CopyFrom(GetX(n)));

			double tol = 1E-13;
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
			DistributedOverlappingIndexer indexer = CreateIndexer(environment);

			Dictionary<int, Vector> localX = environment.CalcNodeData(n => GetX(n));
			var distributedX = new DistributedOverlappingVector(indexer, localX);

			Dictionary<int, Vector> localY = environment.CalcNodeData(n => GetY(n));
			var distributedY = new DistributedOverlappingVector(indexer, localY);

			Dictionary<int, Vector> localZExpected = environment.CalcNodeData(n => 2.0 * GetX(n) + 3.0 * GetY(n));
			var distributedZExpected = new DistributedOverlappingVector(indexer, localZExpected);

			DistributedOverlappingVector distributedZ = distributedX.Copy();
			distributedZ.LinearCombinationIntoThis(2.0, distributedY, 3.0);

			double tol = 1E-13;
			Assert.True(distributedZExpected.Equals(distributedZ, tol));
		}

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

			Dictionary<int, Vector> localX = environment.CalcNodeData(n => GetX(n));
			var distributedX = new DistributedOverlappingVector(indexer, localX);

			Dictionary<int, Vector> localAxExpected = environment.CalcNodeData(n => GetAx(n));
			var distributedAxExpected = new DistributedOverlappingVector(indexer, localAxExpected);

			var distributedAx = new DistributedOverlappingVector(indexer);
			distributedA.MultiplyVector(distributedX, distributedAx);

			double tol = 1E-13;
			Assert.True(distributedAxExpected.Equals(distributedAx, tol));
		}

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

			Dictionary<int, Vector> localAx = environment.CalcNodeData(n => GetAx(n));
			var distributedAx = new DistributedOverlappingVector(indexer, localAx);

			Dictionary<int, Vector> localXExpected = environment.CalcNodeData(n => GetX(n));
			var distributedXExpected = new DistributedOverlappingVector(indexer, localXExpected);

			var pcgBuilder = new PcgAlgorithm.Builder();
			int maxIterations = 12;
			pcgBuilder.MaxIterationsProvider = new FixedMaxIterationsProvider(maxIterations);
			pcgBuilder.ResidualTolerance = 1E-10;
			PcgAlgorithm pcg = pcgBuilder.Build();
			var distributedX = new DistributedOverlappingVector(indexer);
			IterativeStatistics stats = pcg.Solve(distributedA, new IdentityPreconditioner(), distributedAx, distributedX, true);

			double tol = 1E-10;
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

			Dictionary<int, Vector> localAx = environment.CalcNodeData(n => GetAx(n));
			var distributedAx = new DistributedOverlappingVector(indexer, localAx);

			Dictionary<int, Vector> localXExpected = environment.CalcNodeData(n => GetX(n));
			var distributedXExpected = new DistributedOverlappingVector(indexer, localXExpected);

			var pcgBuilder = new BlockPcgAlgorithm.Builder();
			int maxIterations = 12;
			pcgBuilder.MaxIterationsProvider = new FixedMaxIterationsProvider(maxIterations);
			pcgBuilder.ResidualTolerance = 1E-10;
			BlockPcgAlgorithm pcg = pcgBuilder.Build();
			var distributedX = new DistributedOverlappingVector(indexer);
			IterativeStatistics stats = pcg.Solve(distributedA, new IdentityPreconditioner(), distributedAx, distributedX, true);

			double tol = 1E-6;
			Assert.True(stats.HasConverged);
			Assert.True(distributedXExpected.Equals(distributedX, tol));
		}

		[Theory]
		[InlineData(EnvironmentChoice.SequentialSharedEnvironment)]
		[InlineData(EnvironmentChoice.TplSharedEnvironment)]
		public static void TestScaleVectorManaged(EnvironmentChoice env) => TestScaleVector(env.CreateEnvironment());

		internal static void TestScaleVector(IComputeEnvironment environment)
		{
			environment.Initialize(CreateNodeTopology());
			DistributedOverlappingIndexer indexer = CreateIndexer(environment);

			Dictionary<int, Vector> localX = environment.CalcNodeData(n => GetX(n));
			var distributedX = new DistributedOverlappingVector(indexer, localX);

			Dictionary<int, Vector> localZExpected = environment.CalcNodeData(n => -3.0 * GetX(n));
			var distributedZExpected = new DistributedOverlappingVector(indexer, localZExpected);

			DistributedOverlappingVector distributedZ = distributedX.Copy();
			distributedZ.ScaleIntoThis(-3.0);

			double tol = 1E-13;
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
			DistributedOverlappingIndexer indexer = CreateIndexer(environment);

			Dictionary<int, Vector> localInputW = environment.CalcNodeData(n => GetWBeforeSumOverlapping(n));
			var distributedInputW = new DistributedOverlappingVector(indexer, localInputW);

			Dictionary<int, Vector> localOutputW = environment.CalcNodeData(n => GetWAfterSumOverlapping(n));
			var distributedOutputW = new DistributedOverlappingVector(indexer, localOutputW);

			distributedInputW.SumOverlappingEntries();

			double tol = 1E-13;
			Assert.True(distributedOutputW.Equals(distributedInputW, tol));
		}

		private static DistributedOverlappingIndexer CreateIndexer(IComputeEnvironment environment)
		{
			var indexer = new DistributedOverlappingIndexer(environment);
			Action<int> initializeIndexer = n =>
			{
				int previous = n - 1 >= 0 ? n - 1 : numNodes - 1;
				int next = (n + 1) % numNodes;
				var commonEntries = new Dictionary<int, int[]>();
				commonEntries[previous] = new int[] { 0 };
				commonEntries[next] = new int[] { 2 };
				int numEntries = 3;
				indexer.GetLocalComponent(n).Initialize(numEntries, commonEntries);
			};
			environment.DoPerNode(initializeIndexer);
			
			return indexer;
		}

		private static ComputeNodeTopology CreateNodeTopology() //TODOMPI: this is also described in Environments.Tests. 
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

		private static Vector GetAx(int nodeID)
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

		private static Matrix GetMatrixA(int nodeID)
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

		private static Vector GetX(int nodeID)
		{
			var x = Vector.CreateZero(3);
			for (int i = 0; i < x.Length; ++i)
			{
				x[i] = 2 * nodeID + i;
			}
			if (nodeID == numNodes - 1)
			{
				x[x.Length - 1] = 0;
			}
			return x;
		}

		private static Vector GetY(int nodeID)
		{
			Vector y = GetX(nodeID);
			y.DoToAllEntriesIntoThis(a => 10 + a);
			return y;
		}

		private static double GetXDotY() => 1166;

		private static Vector GetWBeforeSumOverlapping(int nodeID)
		{
			var w = Vector.CreateZero(3);
			for (int i = 0; i < w.Length; ++i)
			{
				w[i] = 3 * nodeID + i;
			}
			return w;
		}

		private static Vector GetWAfterSumOverlapping(int nodeID)
		{
			int previous = nodeID - 1 >= 0 ? nodeID - 1 : numNodes - 1;
			int next = (nodeID + 1) % numNodes;

			Vector w = GetWBeforeSumOverlapping(nodeID);
			Vector wPrevious = GetWBeforeSumOverlapping(previous);
			Vector wNext = GetWBeforeSumOverlapping(next);

			w[0] += wPrevious[2];
			w[2] += wNext[0];

			return w;
		}
	}
}
