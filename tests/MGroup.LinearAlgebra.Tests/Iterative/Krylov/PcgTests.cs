namespace MGroup.LinearAlgebra.Tests.Iterative.Krylov
{
	using MGroup.LinearAlgebra.Implementations;
	using MGroup.LinearAlgebra.Iterative;
	using MGroup.LinearAlgebra.Iterative.ConjugateGradient;
	using MGroup.LinearAlgebra.Iterative.PreconditionedConjugateGradient;
	using MGroup.LinearAlgebra.Iterative.Preconditioning;
	using MGroup.LinearAlgebra.Iterative.Termination.Iterations;
	using MGroup.LinearAlgebra.Matrices;
	using MGroup.LinearAlgebra.Tests.TestData;
	using MGroup.LinearAlgebra.Tests.Utilities;
	using MGroup.LinearAlgebra.Vectors;

	using Xunit;

	/// <summary>
	/// Tests for <see cref="PcgAlgorithm"/>.
	/// </summary>
	public static class PcgTests
	{
		private static readonly MatrixComparer comparer = new MatrixComparer(1E-5);

		[Theory]
		[MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
		private static void TestPosDefDenseSystem(IImplementationProvider provider)
		{
			TestSettings.RunMultiproviderTest(provider, delegate ()
			{
				var A = Matrix.CreateFromArray(SymmPosDef10by10.Matrix);
				var b = Vector.CreateFromArray(SymmPosDef10by10.Rhs);
				var xExpected = Vector.CreateFromArray(SymmPosDef10by10.Lhs);

				var builder = new PcgAlgorithm.Factory();
				builder.ResidualTolerance = 1E-7;
				builder.MaxIterationsProvider = new PercentageMaxIterationsProvider(1.0);
				var pcg = builder.Build();
				var M = new JacobiPreconditioner();
				M.UpdateMatrix(A, true);
				var xComputed = Vector.CreateZero(A.NumRows);
				var stats = pcg.Solve(A, M, b, xComputed, true, () => Vector.CreateZero(b.Length));
				comparer.AssertEqual(xExpected, xComputed);
			});
		}

		[Theory]
		[MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
		private static void TestPosDefSparseSystem(IImplementationProvider provider)
		{
			TestSettings.RunMultiproviderTest(provider, delegate ()
			{
				var A = Matrix.CreateFromArray(SparsePosDef10by10.Matrix);
				var b = Vector.CreateFromArray(SparsePosDef10by10.Rhs);
				var xExpected = Vector.CreateFromArray(SparsePosDef10by10.Lhs);

				var builder = new PcgAlgorithm.Factory();
				builder.ResidualTolerance = 1E-7;
				builder.MaxIterationsProvider = new PercentageMaxIterationsProvider(1.0);
				var pcg = builder.Build();
				var M = new JacobiPreconditioner();
				M.UpdateMatrix(A, true);
				var xComputed = Vector.CreateZero(A.NumRows);
				var stats = pcg.Solve(A, M, b, xComputed, true, () => Vector.CreateZero(b.Length));
				comparer.AssertEqual(xExpected, xComputed);
			});
		}

		[Theory]
		[MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
		private static void TestIndefiniteSystem(IImplementationProvider provider)
		{
			TestSettings.RunMultiproviderTest(provider, delegate ()
			{
				(var A, var b, var xExpected, var M) = DiagonalIndefinite.BuildIndefiniteSystem(20);
				var builder = new CGAlgorithm.Builder();
				builder.ResidualTolerance = 1E-6;
				builder.MaxIterationsProvider = new PercentageMaxIterationsProvider(1.0);
				var cg = builder.Build();
				var xComputed = Vector.CreateZero(A.NumRows);
				var stats = cg.Solve(A, b, xComputed, true);
				Assert.False(comparer.AreEqual(xExpected, xComputed));
			});
		}
	}
}
