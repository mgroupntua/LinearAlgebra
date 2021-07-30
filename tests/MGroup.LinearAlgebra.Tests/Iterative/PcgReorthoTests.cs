using System.Diagnostics;
using MGroup.LinearAlgebra.Iterative;
using MGroup.LinearAlgebra.Iterative.PreconditionedConjugateGradient;
using MGroup.LinearAlgebra.Iterative.Preconditioning;
using MGroup.LinearAlgebra.Iterative.Termination;
using MGroup.LinearAlgebra.Matrices;
using MGroup.LinearAlgebra.Tests.TestData;
using MGroup.LinearAlgebra.Tests.Utilities;
using MGroup.LinearAlgebra.Vectors;
using Xunit;

namespace MGroup.LinearAlgebra.Tests.Iterative
{
	/// <summary>
	/// Tests for <see cref="ReorthogonalizedPcg"/>.
	/// Authors: Serafeim Bakalakos
	/// </summary>
	public static class PcgReorthoTests
	{
		private static readonly MatrixComparer comparer = new MatrixComparer(1E-5);

		[Theory]
		[InlineData(0.1, 5, 10)]
		[InlineData(0.01, 5, 20)]
		private static void TestNearbyProblems(double noiseWidth, int maxIterations, int numRhsVectors)
		{
			int order = SymmPosDef10by10.Order;
			var A = Matrix.CreateFromArray(SymmPosDef10by10.Matrix);
			var builder = new ReorthogonalizedPcg.Builder();
			builder.ResidualTolerance = 1E-6;
			builder.MaxIterationsProvider = new PercentageMaxIterationsProvider(1.0);
			builder.Convergence = new RhsNormalizedConvergence();
			var pcg = builder.Build();
			var M = new JacobiPreconditioner(A.GetDiagonalAsArray());

			// Initial run
			Vector x0 = Vector.CreateWithValue(order, 1);
			Vector x0Expected = x0.Copy();
			Vector b0 = A * x0Expected;
			Vector x0Computed = Vector.CreateZero(A.NumRows);
			IterativeStatistics stats0 = pcg.Solve(A, M, b0, x0Computed, true, () => Vector.CreateZero(order));
			 Debug.WriteLine($"Initial run: iterations = {stats0.NumIterationsRequired}");
			comparer.AssertEqual(x0Expected, x0Computed);

			// Subsequent runs
			int seed = 345;
			for (int i = 0; i < numRhsVectors; ++i)
			{
				Vector dx = Vector.CreateFromArray(RandomMatrices.CreateRandomVector(order, seed));
				Vector xExpected = x0 + noiseWidth * dx;
				Vector b = A * xExpected;

				pcg.Clear(); //TODO: preferably do not call this.
				//pcg.ReorthoCache.Clear();

				Vector xComputed = Vector.CreateZero(A.NumRows);
				IterativeStatistics stats = pcg.Solve(A, M, b, xComputed, true, () => Vector.CreateZero(b.Length));
				Debug.WriteLine($"Subsequent run: iterations = {stats.NumIterationsRequired}");
				comparer.AssertEqual(xExpected, xComputed);
				Assert.InRange(stats.NumIterationsRequired, 1, maxIterations);
			}
		}

		[Fact]
		private static void TestPosDefDenseSystem()
		{
			var A = Matrix.CreateFromArray(SymmPosDef10by10.Matrix);
			var b = Vector.CreateFromArray(SymmPosDef10by10.Rhs);
			var xExpected = Vector.CreateFromArray(SymmPosDef10by10.Lhs);

			var builder = new ReorthogonalizedPcg.Builder();
			builder.ResidualTolerance = 1E-7;
			builder.MaxIterationsProvider = new PercentageMaxIterationsProvider(1.0);
			var pcg = builder.Build();
			var M = new JacobiPreconditioner(A.GetDiagonalAsArray());
			Vector xComputed = Vector.CreateZero(A.NumRows);
			IterativeStatistics stats = pcg.Solve(A, M, b, xComputed, true, () => Vector.CreateZero(b.Length));
			comparer.AssertEqual(xExpected, xComputed);
		}

		[Fact]
		private static void TestPosDefSparseSystem()
		{
			var A = Matrix.CreateFromArray(SparsePosDef10by10.Matrix);
			var b = Vector.CreateFromArray(SparsePosDef10by10.Rhs);
			var xExpected = Vector.CreateFromArray(SparsePosDef10by10.Lhs);

			var builder = new ReorthogonalizedPcg.Builder();
			builder.ResidualTolerance = 1E-7;
			builder.MaxIterationsProvider = new PercentageMaxIterationsProvider(1.0);
			var pcg = builder.Build();
			var M = new JacobiPreconditioner(A.GetDiagonalAsArray());
			Vector xComputed = Vector.CreateZero(A.NumRows);
			IterativeStatistics stats = pcg.Solve(A, M, b, xComputed, true, () => Vector.CreateZero(b.Length));
			comparer.AssertEqual(xExpected, xComputed);
		}
	}
}
