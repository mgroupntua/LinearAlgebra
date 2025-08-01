using MGroup.LinearAlgebra.Implementations;
using MGroup.LinearAlgebra.Iterative.MinimumResidual;
using MGroup.LinearAlgebra.Iterative.Preconditioning;
using MGroup.LinearAlgebra.Matrices;
using MGroup.LinearAlgebra.Tests.TestData;
using MGroup.LinearAlgebra.Tests.Utilities;
using MGroup.LinearAlgebra.Vectors;

using Xunit;

namespace MGroup.LinearAlgebra.Tests.Iterative.Krylov
{
	/// <summary>
	/// Tests for <see cref="SparsityPatternSymmetric"/>.
	/// </summary>
	public static class MinResTests
	{
		[Theory]
		[MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
		public static void TestIndefiniteSystem(IImplementationProvider provider)
		{
			TestSettings.RunMultiproviderTest(provider, delegate ()
			{
				var comparer = new MatrixComparer(1E-4);
				var residualTolerance = 1e-8;
				(var A, var b, var xExpected, var M) = DiagonalIndefinite.BuildIndefiniteSystem(20);
				var minres = new MinRes(A.NumRows, residualTolerance, 0, false, false);
				(var xComputed, var stats) = minres.Solve(A, b);
				comparer.AssertEqual(xExpected, xComputed);
			});
		}

		[Theory]
		[MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
		public static void TestPosDefDenseSystem(IImplementationProvider provider)
		{
			TestSettings.RunMultiproviderTest(provider, delegate ()
			{
				var comparer = new MatrixComparer(1E-6);
				var n = SparsePosDef10by10.Order;
				var A = Matrix.CreateFromArray(SparsePosDef10by10.Matrix);
				var b = Vector.CreateFromArray(SparsePosDef10by10.Rhs);
				var xExpected = Vector.CreateFromArray(SparsePosDef10by10.Lhs);
				var minres = new MinRes(n, 1e-10, 0, false, false);
				(var xComputed, var stats) = minres.Solve(A, b);
				comparer.AssertEqual(xExpected, xComputed);
			});
		}

		[Theory]
		[MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
		public static void TestPosDefSparseSystem(IImplementationProvider provider)
		{
			TestSettings.RunMultiproviderTest(provider, delegate ()
			{
				var comparer = new MatrixComparer(1E-6);
				var n = SparsePosDef10by10.Order;
				var A = Matrix.CreateFromArray(SparsePosDef10by10.Matrix);
				var b = Vector.CreateFromArray(SparsePosDef10by10.Rhs);
				var xExpected = Vector.CreateFromArray(SparsePosDef10by10.Lhs);
				var minres = new MinRes(n, 1e-10, 0, false, false);
				(var xComputed, var stats) = minres.Solve(A, b);
				comparer.AssertEqual(xExpected, xComputed);
			});
		}
	}
}
