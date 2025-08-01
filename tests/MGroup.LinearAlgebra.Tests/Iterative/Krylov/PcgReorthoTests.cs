using System.Collections.Generic;
using System.Diagnostics;

using MGroup.LinearAlgebra.Iterative;
using MGroup.LinearAlgebra.Iterative.PreconditionedConjugateGradient;
using MGroup.LinearAlgebra.Iterative.PreconditionedConjugateGradient.Reorthogonalization;
using MGroup.LinearAlgebra.Iterative.Preconditioning;
using MGroup.LinearAlgebra.Iterative.Termination.Iterations;
using MGroup.LinearAlgebra.Matrices;
using MGroup.LinearAlgebra.Tests.TestData;
using MGroup.LinearAlgebra.Tests.TestData.FiniteElementMatrices;
using MGroup.LinearAlgebra.Tests.Utilities;
using MGroup.LinearAlgebra.Vectors;

using Xunit;

namespace MGroup.LinearAlgebra.Tests.Iterative.Krylov
{
	/// <summary>
	/// Tests for <see cref="ReorthogonalizedPcg"/>.
	/// </summary>
	public static class PcgReorthoTests
	{
		private static readonly MatrixComparer comparer = new MatrixComparer(1E-5);

		//[Fact]
		public static void InvestigateNoiseStagnation()
		{
			double noiseWidth = 100;

			var order = 100;
			//var A = Matrix.CreateFromArray(MultiDiagonalMatrices.CreateSymmetric(order, new int[] { 0, 1, 5, 7, 12 }));
			var valueOfEachDiagonal = new Dictionary<int, double>();
			valueOfEachDiagonal[0] = 10.0;
			valueOfEachDiagonal[1] = 4.0;
			valueOfEachDiagonal[5] = 3.0;
			var A = Matrix.CreateFromArray(MultiDiagonalMatrices.CreateSymmetric(order, valueOfEachDiagonal));
			var M = new IdentityPreconditioner();
			//var M = new JacobiPreconditioner(A.GetDiagonalAsArray());

			// Method A: Regular PCG
			var pcgBuilder = new PcgAlgorithm.Factory();
			pcgBuilder.ResidualTolerance = 1E-15;
			pcgBuilder.MaxIterationsProvider = new FixedMaxIterationsProvider(50);
			var pcg = pcgBuilder.Build();

			// Method B: Reorthogonalized PCG, but without keeping direction vectors from previous solutions.
			var pcgReorthoRestartBuilder = new ReorthogonalizedPcg.Factory();
			pcgReorthoRestartBuilder.ResidualTolerance = 1E-15;
			pcgReorthoRestartBuilder.MaxIterationsProvider = new FixedMaxIterationsProvider(50);
			var pcgReorthoRestart = pcgReorthoRestartBuilder.Build();

			// Method C: Reorthogonalized PCG, where the second solution will reuse direction vectors from the first
			var pcgReorthoBuilder = new ReorthogonalizedPcg.Factory();
			pcgReorthoBuilder.ResidualTolerance = 1E-15;
			pcgReorthoBuilder.MaxIterationsProvider = new FixedMaxIterationsProvider(50);
			var pcgReortho = pcgReorthoBuilder.Build();

			// Initial rhs
			var x0 = Vector.CreateWithValue(order, 1);
			var x0Expected = x0.Copy();
			var b0 = A * x0Expected;

			var xA = Vector.CreateZero(A.NumRows);
			var statsA = pcg.Solve(A, M, b0, xA, true);
			Debug.WriteLine($"Initial run - method A: iterations = {statsA.NumIterationsRequired}");

			var xB = Vector.CreateZero(A.NumRows);
			var statsB = pcgReorthoRestart.Solve(A, M, b0, xB, true);
			Debug.WriteLine($"Initial run - method B iterations = {statsB.NumIterationsRequired}");

			var xC = Vector.CreateZero(A.NumRows);
			var statsC = pcgReortho.Solve(A, M, b0, xC, true);
			Debug.WriteLine($"Initial run - method C: iterations = {statsC.NumIterationsRequired}");

			// Perturbed rhs
			var seed = 345;
			var dx = Vector.CreateFromArray(RandomMatrices.CreateRandomVector(order, seed));

			var x1Expected = x0 + noiseWidth * dx;
			var b1 = A * x1Expected;

			xA = Vector.CreateZero(A.NumRows);
			statsA = pcg.Solve(A, M, b1, xA, true);
			Debug.WriteLine($"2nd run, noise = {noiseWidth} - method A: iterations = {statsA.NumIterationsRequired}");

			xB = Vector.CreateZero(A.NumRows);
			pcgReorthoRestart.ReorthoCache.Clear();
			statsB = pcgReorthoRestart.Solve(A, M, b1, xB, true);
			Debug.WriteLine($"2nd run, noise = {noiseWidth} - method B iterations = {statsB.NumIterationsRequired}");

			xC = Vector.CreateZero(A.NumRows);
			statsC = pcgReortho.Solve(A, M, b1, xC, true);
			Debug.WriteLine($"2nd run, noise = {noiseWidth} - method C: iterations = {statsC.NumIterationsRequired}");
		}

		//[Fact]
		public static void InvestigatePFetiDPCoarseProblem2D()
		{
			var order = PFetiDPCoarseProblem2D.Order;
			var A = Matrix.CreateFromArray(PFetiDPCoarseProblem2D.MatrixScc);
			var M = new IdentityPreconditioner();

			// Method A: Regular PCG
			var pcgBuilder = new PcgAlgorithm.Factory();
			pcgBuilder.ResidualTolerance = 1E-20;
			pcgBuilder.MaxIterationsProvider = new FixedMaxIterationsProvider(50);
			var pcg = pcgBuilder.Build();

			// Method B: Reorthogonalized PCG, but without keeping direction vectors from previous solutions.
			var pcgReorthoRestartBuilder = new ReorthogonalizedPcg.Factory();
			pcgReorthoRestartBuilder.ResidualTolerance = 1E-20;
			pcgReorthoRestartBuilder.MaxIterationsProvider = new FixedMaxIterationsProvider(50);
			var pcgReorthoRestart = pcgReorthoRestartBuilder.Build();

			// Method C: Reorthogonalized PCG, where the second solution will reuse direction vectors from the first
			var pcgReorthoBuilder = new ReorthogonalizedPcg.Factory();
			pcgReorthoBuilder.ResidualTolerance = 1E-20;
			pcgReorthoBuilder.MaxIterationsProvider = new FixedMaxIterationsProvider(50);
			var pcgReortho = pcgReorthoBuilder.Build();

			// Initial rhs
			var b = Vector.CreateFromArray(PFetiDPCoarseProblem2D.RhsVectors[0]);
			var xExpected = Vector.CreateFromArray(PFetiDPCoarseProblem2D.SolutionVectors[0]);

			var xA = Vector.CreateZero(A.NumRows);
			var statsA = pcg.Solve(A, M, b, xA, true);
			Assert.True(xExpected.Equals(xA, 1E-10));
			Debug.WriteLine($"Initial run - method A: iterations = {statsA.NumIterationsRequired}");

			var xB = Vector.CreateZero(A.NumRows);
			var statsB = pcgReorthoRestart.Solve(A, M, b, xB, true);
			Assert.True(xExpected.Equals(xB, 1E-10));
			Debug.WriteLine($"Initial run - method B iterations = {statsB.NumIterationsRequired}");

			var xC = Vector.CreateZero(A.NumRows);
			var statsC = pcgReortho.Solve(A, M, b, xC, true);
			Assert.True(xExpected.Equals(xC, 1E-10));
			Debug.WriteLine($"Initial run - method C: iterations = {statsC.NumIterationsRequired}");

			// Next rhs
			b = Vector.CreateFromArray(PFetiDPCoarseProblem2D.RhsVectors[1]);
			xExpected = Vector.CreateFromArray(PFetiDPCoarseProblem2D.SolutionVectors[1]);

			xA = Vector.CreateZero(A.NumRows);
			statsA = pcg.Solve(A, M, b, xA, true);
			Assert.True(xExpected.Equals(xA, 1E-10));
			Debug.WriteLine($"Initial run - method A: iterations = {statsA.NumIterationsRequired}");

			xB = Vector.CreateZero(A.NumRows);
			pcgReorthoRestart.ReorthoCache.Clear();
			statsB = pcgReorthoRestart.Solve(A, M, b, xB, true);
			Assert.True(xExpected.Equals(xB, 1E-10));
			Debug.WriteLine($"Initial run - method B iterations = {statsB.NumIterationsRequired}");

			xC = Vector.CreateZero(A.NumRows);
			statsC = pcgReortho.Solve(A, M, b, xC, true);
			Assert.True(xExpected.Equals(xC, 1E-10));
			Debug.WriteLine($"Initial run - method C: iterations = {statsC.NumIterationsRequired}");
		}

		[Theory]
		[InlineData(0.1, 5, 10)]
		[InlineData(0.01, 5, 20)]
		public static void TestNearbyProblems(double noiseWidth, int maxIterations, int numRhsVectors)
		{
			var order = SymmPosDef10by10.Order;
			var A = Matrix.CreateFromArray(SymmPosDef10by10.Matrix);
			var builder = new ReorthogonalizedPcg.Factory();
			builder.ResidualTolerance = 1E-6;
			builder.MaxIterationsProvider = new PercentageMaxIterationsProvider(1.0);
			builder.Convergence = new RhsNormalizedConvergence();
			var pcg = builder.Build();
			var M = new JacobiPreconditioner();
			M.UpdateMatrix(A, true);

			// Initial run
			var x0 = Vector.CreateWithValue(order, 1);
			var x0Expected = x0.Copy();
			var b0 = A * x0Expected;
			var x0Computed = Vector.CreateZero(A.NumRows);
			var stats0 = pcg.Solve(A, M, b0, x0Computed, true);
			Debug.WriteLine($"Initial run: iterations = {stats0.NumIterationsRequired}");
			comparer.AssertEqual(x0Expected, x0Computed);

			// Subsequent runs
			var seed = 345;
			for (var i = 0; i < numRhsVectors; ++i)
			{
				var dx = Vector.CreateFromArray(RandomMatrices.CreateRandomVector(order, seed));
				var xExpected = x0 + noiseWidth * dx;
				var b = A * xExpected;

				pcg.Clear(); //TODO: preferably do not call this.
							 //pcg.ReorthoCache.Clear();

				var xComputed = Vector.CreateZero(A.NumRows);
				var stats = pcg.Solve(A, M, b, xComputed, true);
				Debug.WriteLine($"Subsequent run: iterations = {stats.NumIterationsRequired}");
				comparer.AssertEqual(xExpected, xComputed);
				Assert.InRange(stats.NumIterationsRequired, 1, maxIterations);
			}
		}

		[Fact]
		public static void TestPosDefDenseSystem()
		{
			var A = Matrix.CreateFromArray(SymmPosDef10by10.Matrix);
			var b = Vector.CreateFromArray(SymmPosDef10by10.Rhs);
			var xExpected = Vector.CreateFromArray(SymmPosDef10by10.Lhs);

			var factory = new ReorthogonalizedPcg.Factory();
			factory.ResidualTolerance = 1E-7;
			factory.MaxIterationsProvider = new PercentageMaxIterationsProvider(1.0);
			var pcg = factory.Build();
			var M = new JacobiPreconditioner();
			M.UpdateMatrix(A, true);
			var xComputed = Vector.CreateZero(A.NumRows);
			var stats = pcg.Solve(A, M, b, xComputed, true);
			comparer.AssertEqual(xExpected, xComputed);
		}

		[Fact]
		public static void TestPosDefSparseSystem()
		{
			var A = Matrix.CreateFromArray(SparsePosDef10by10.Matrix);
			var b = Vector.CreateFromArray(SparsePosDef10by10.Rhs);
			var xExpected = Vector.CreateFromArray(SparsePosDef10by10.Lhs);

			var factory = new ReorthogonalizedPcg.Factory();
			factory.ResidualTolerance = 1E-7;
			factory.MaxIterationsProvider = new PercentageMaxIterationsProvider(1.0);
			var pcg = factory.Build();
			var M = new JacobiPreconditioner();
			M.UpdateMatrix(A, true);
			var xComputed = Vector.CreateZero(A.NumRows);
			var stats = pcg.Solve(A, M, b, xComputed, true);
			comparer.AssertEqual(xExpected, xComputed);
		}
	}
}
