namespace MGroup.LinearAlgebra.Tests.Iterative.Krylov
{
	using System.Diagnostics;

	using MGroup.LinearAlgebra.Exceptions;
	using MGroup.LinearAlgebra.Implementations;
	using MGroup.LinearAlgebra.Iterative;
	using MGroup.LinearAlgebra.Iterative.ConjugateGradient;
	using MGroup.LinearAlgebra.Iterative.PreconditionedConjugateGradient;
	using MGroup.LinearAlgebra.Iterative.PreconditionedConjugateGradient.Logging;
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
		public static void TestPosDefDenseSystem(IImplementationProvider provider)
		{
			TestSettings.RunMultiproviderTest(provider, delegate ()
			{
				var A = Matrix.CreateFromArray(SymmPosDef10by10.Matrix);
				var b = Vector.CreateFromArray(SymmPosDef10by10.Rhs);
				var xExpected = Vector.CreateFromArray(SymmPosDef10by10.Lhs);

				var factory = new PcgAlgorithm.Factory();
				factory.ResidualTolerance = 1E-7;
				factory.MaxIterationsProvider = new PercentageMaxIterationsProvider(1.0);
				var pcg = factory.Build();
				var M = new JacobiPreconditioner();
				M.UpdateMatrix(A, true);

				var xComputed = Vector.CreateZero(A.NumRows);
				IterativeStatistics stats = pcg.Solve(A, M, b, xComputed, true);
				comparer.AssertEqual(xExpected, xComputed);
			});
		}

		[Theory]
		[MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
		public static void TestPosDefSparseSystem(IImplementationProvider provider)
		{
			TestSettings.RunMultiproviderTest(provider, delegate ()
			{
				var A = Matrix.CreateFromArray(SparsePosDef10by10.Matrix);
				var b = Vector.CreateFromArray(SparsePosDef10by10.Rhs);
				var xExpected = Vector.CreateFromArray(SparsePosDef10by10.Lhs);

				var factory = new PcgAlgorithm.Factory();
				factory.ResidualTolerance = 1E-7;
				factory.MaxIterationsProvider = new PercentageMaxIterationsProvider(1.0);
				factory.Logger = new ResidualNormRatioLogger();
				var pcg = factory.Build();
				var M = new JacobiPreconditioner();
				M.UpdateMatrix(A, true);

				var xComputed = Vector.CreateZero(A.NumRows);
				IterativeStatistics stats = pcg.Solve(A, M, b, xComputed, true);
				comparer.AssertEqual(xExpected, xComputed);
				Debug.WriteLine(pcg.Logger.Report());
			});
		}

		[Theory]
		[MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
		public static void TestIndefiniteSystem(IImplementationProvider provider)
		{
			TestSettings.RunMultiproviderTest(provider, delegate ()
			{
				(Matrix A, Vector b, Vector xExpected, IPreconditioner M) = DiagonalIndefinite.BuildIndefiniteSystem(20);
				var factory = new PcgAlgorithm.Factory();
				factory.ResidualTolerance = 1E-6;
				factory.MaxIterationsProvider = new PercentageMaxIterationsProvider(1.0);
				factory.ThrowExceptionIfNotConvergence = true;
				var pcg = factory.Build();

				var xComputed = Vector.CreateZero(A.NumRows);
				Assert.Throws<IterativeMethodDidNotConvergeException>(() =>
				{
					IterativeStatistics stats = pcg.Solve(A, M, b, xComputed, true);
				});
			});
		}
	}
}
