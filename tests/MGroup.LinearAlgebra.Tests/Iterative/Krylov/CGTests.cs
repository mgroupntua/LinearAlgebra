using MGroup.LinearAlgebra.Implementations;
using MGroup.LinearAlgebra.Iterative;
using MGroup.LinearAlgebra.Iterative.ConjugateGradient;
using MGroup.LinearAlgebra.Iterative.Termination.Iterations;
using MGroup.LinearAlgebra.Matrices;
using MGroup.LinearAlgebra.Tests.TestData;
using MGroup.LinearAlgebra.Tests.Utilities;
using MGroup.LinearAlgebra.Vectors;

using Xunit;

namespace MGroup.LinearAlgebra.Tests.Iterative.Krylov
{
	/// <summary>
	/// Tests for <see cref="PcgAlgorithm"/>.
	/// </summary>
	public static class CGTests
	{
		private static readonly MatrixComparer comparer = new MatrixComparer(1E-5);

		[Theory]
		[MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
		private static void TestDenseSystem(IImplementationProvider provider)
		{
			TestSettings.RunMultiproviderTest(provider, delegate ()
			{
				var A = Matrix.CreateFromArray(SymmPosDef10by10.Matrix);
				var b = Vector.CreateFromArray(SymmPosDef10by10.Rhs);
				var xExpected = Vector.CreateFromArray(SymmPosDef10by10.Lhs);

				var builder = new CGAlgorithm.Builder();
				builder.ResidualTolerance = 1E-7;
				builder.MaxIterationsProvider = new PercentageMaxIterationsProvider(1.0);
				var cg = builder.Build();
				var xComputed = Vector.CreateZero(A.NumRows);
				var stats = cg.Solve(A, b, xComputed, true);
				comparer.AssertEqual(xExpected, xComputed);
			});
		}

		[Theory]
		[MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
		private static void TestSparseSystem(IImplementationProvider provider)
		{
			TestSettings.RunMultiproviderTest(provider, delegate ()
			{
				var A = Matrix.CreateFromArray(SparsePosDef10by10.Matrix);
				var b = Vector.CreateFromArray(SparsePosDef10by10.Rhs);
				var xExpected = Vector.CreateFromArray(SparsePosDef10by10.Lhs);

				var builder = new CGAlgorithm.Builder();
				builder.ResidualTolerance = 1E-7;
				builder.MaxIterationsProvider = new PercentageMaxIterationsProvider(1.0);
				var cg = builder.Build();
				var xComputed = Vector.CreateZero(A.NumRows);
				var stats = cg.Solve(A, b, xComputed, true);
				comparer.AssertEqual(xExpected, xComputed);
			});
		}
	}
}
