using System.Text;

using MGroup.LinearAlgebra.Iterative;
using MGroup.LinearAlgebra.Iterative.ConjugateGradient;
using MGroup.LinearAlgebra.Iterative.ConjugateGradient;
using MGroup.LinearAlgebra.Iterative.Termination;
using MGroup.LinearAlgebra.Matrices;
using MGroup.LinearAlgebra.Output;
using MGroup.LinearAlgebra.Tests.TestData;
using MGroup.LinearAlgebra.Tests.Utilities;
using MGroup.LinearAlgebra.Vectors;
using Xunit;

namespace MGroup.LinearAlgebra.Tests.Iterative
{
	/// <summary>
	/// Tests for <see cref="GaussSeidel"/>.
	/// Authors: Gerasimos Sotiropoulos 
	/// </summary>
	public static class GaussSeidelTests
	{
		private static readonly MatrixComparer comparer = new MatrixComparer(1E-5);


		[Theory]
		[InlineData(true)]
		[InlineData(false)]
		private static void TestSparseSystem(bool forwardGaussSeidel)
		{
			var A = CsrMatrix.CreateFromArrays(SparsePosDef10by10.Order, SparsePosDef10by10.Order, SparsePosDef10by10.CsrValues, SparsePosDef10by10.CsrColIndices, SparsePosDef10by10.CsrRowOffsets, true);
			var b = Vector.CreateFromArray(SparsePosDef10by10.Rhs);
			var xExpected = Vector.CreateFromArray(SparsePosDef10by10.Lhs);

			var builder = new GaussSeidel.Builder();
			builder.ConvergenceTolerance = 1E-7;
			builder.MaxIterationsProvider = new FixedMaxIterationsProvider(100);
			builder.ForwardGaussSeidel = true;
			var gs = builder.Build();
			var xComputed = Vector.CreateZero(A.NumRows);

			IterativeStatistics stats = gs.Solve(A, b, xComputed);
			comparer.AssertEqual(xExpected, xComputed);
		}
	}
}
