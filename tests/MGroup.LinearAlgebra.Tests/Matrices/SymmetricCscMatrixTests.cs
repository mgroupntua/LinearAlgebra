using MGroup.LinearAlgebra.Commons;
using MGroup.LinearAlgebra.Matrices;
using MGroup.LinearAlgebra.Tests.TestData;
using MGroup.LinearAlgebra.Tests.Utilities;
using MGroup.LinearAlgebra.Vectors;

using Xunit;

namespace MGroup.LinearAlgebra.Tests.Matrices
{
	/// <summary>
	/// Tests for <see cref="SymmetricCscMatrix"/>.
	/// Authors: Serafeim Bakalakos
	/// </summary>
	public static class SymmetricCscMatrixTests
	{
		private static readonly MatrixComparer comparer = new MatrixComparer(1E-13);

		[Fact]
		private static void TestEquality()
		{
			var full = Matrix.CreateFromArray(SparsePosDef10by10.Matrix);
			var cscSymm = SymmetricCscMatrix.CreateFromArrays(
				SparsePosDef10by10.Order, SparsePosDef10by10.SymmetricCscValues, 
				SparsePosDef10by10.SymmetricCscRowIndices, SparsePosDef10by10.SymmetricCscColOffsets,
				true);
			Assert.True(cscSymm.Equals(full));
		}

		[Theory]
		[MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
		private static void TestMatrixVectorMultiplication(LinearAlgebraProviderChoice providers)
		{
			TestSettings.RunMultiproviderTest(providers, delegate ()
			{
				var A = SymmetricCscMatrix.CreateFromArrays(
					SparsePosDef10by10.Order, SparsePosDef10by10.SymmetricCscValues,
					SparsePosDef10by10.SymmetricCscRowIndices, SparsePosDef10by10.SymmetricCscColOffsets,
					true);

				// Multiply() - untransposed
				var x = Vector.CreateFromArray(SparsePosDef10by10.Lhs);
				var bExpected = Vector.CreateFromArray(SparsePosDef10by10.Rhs);
				IVector bComputed = A.Multiply(x, false);
				comparer.AssertEqual(bExpected, bComputed);

				// Multiply() - transposed
				x = Vector.CreateFromArray(SparsePosDef10by10.Lhs);
				bExpected = Vector.CreateFromArray(SparsePosDef10by10.Rhs); // same as before, since the matrix is symmetric
				bComputed = A.Multiply(x, true);
				comparer.AssertEqual(bExpected, bComputed);

				// MultiplyRight() - untransposed
				x = Vector.CreateFromArray(SparsePosDef10by10.Lhs);
				bExpected = Vector.CreateFromArray(SparsePosDef10by10.Rhs);
				Vector bDenseComputed = A.MultiplyRight(x, false);
				comparer.AssertEqual(bExpected, bDenseComputed);

				// MultiplyRight() - transposed
				x = Vector.CreateFromArray(SparsePosDef10by10.Lhs);
				bExpected = Vector.CreateFromArray(SparsePosDef10by10.Rhs); // same as before, since the matrix is symmetric
				bDenseComputed = A.MultiplyRight(x, true);
				comparer.AssertEqual(bExpected, bDenseComputed);
			});
		}

		[Theory]
		[MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
		private static void TestMatrixVectorMultiplicationIntoResult(LinearAlgebraProviderChoice providers)
		{
			TestSettings.RunMultiproviderTest(providers, delegate ()
			{
				// The result vectors will first be set to some non zero values to make sure that the result overwrites 
				// them instead of being added to them.
				var A = SymmetricCscMatrix.CreateFromArrays(
					SparsePosDef10by10.Order, SparsePosDef10by10.SymmetricCscValues,
					SparsePosDef10by10.SymmetricCscRowIndices, SparsePosDef10by10.SymmetricCscColOffsets,
					true);

				// MultiplyIntoResult() - untransposed
				var x = Vector.CreateFromArray(SparsePosDef10by10.Lhs);
				var bExpected = Vector.CreateFromArray(SparsePosDef10by10.Rhs);
				Vector bComputed = Vector.CreateWithValue(SparsePosDef10by10.Order, 1.0);
				A.MultiplyIntoResult(x, bComputed, false);
				comparer.AssertEqual(bExpected, bComputed);

				// MultiplyIntoResult() - transposed
				x = Vector.CreateFromArray(SparsePosDef10by10.Lhs);
				bExpected = Vector.CreateFromArray(SparsePosDef10by10.Rhs);
				bComputed = Vector.CreateWithValue(SparsePosDef10by10.Order, 1.0);
				A.MultiplyIntoResult(x, bComputed, true);
				comparer.AssertEqual(bExpected, bComputed);
			});
		}
	}
}
