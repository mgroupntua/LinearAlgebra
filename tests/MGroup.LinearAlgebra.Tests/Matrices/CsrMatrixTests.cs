using MGroup.LinearAlgebra.Commons;
using MGroup.LinearAlgebra.Exceptions;
using MGroup.LinearAlgebra.Matrices;
using MGroup.LinearAlgebra.Tests.TestData;
using MGroup.LinearAlgebra.Tests.Utilities;
using MGroup.LinearAlgebra.Vectors;
using Xunit;

namespace MGroup.LinearAlgebra.Tests.Matrices
{
    /// <summary>
    /// Tests for <see cref="CsrMatrix"/>.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class CsrMatrixTests
    {
        private static readonly MatrixComparer comparer = new MatrixComparer(1E-13);

        [Fact]
        private static void TestClear()
        {
            var zero = Matrix.CreateZero(SparseRectangular10by5.NumRows, SparseRectangular10by5.NumCols);
            var csr = CsrMatrix.CreateFromArrays(SparseRectangular10by5.NumRows, SparseRectangular10by5.NumCols,
                SparseRectangular10by5.CsrValues, SparseRectangular10by5.CsrColIndices, SparseRectangular10by5.CsrRowOffsets,
                true);
            csr.Clear();
            comparer.AssertEqual(zero, csr);
        }

		[Fact]
		private static void TestCreateFromDense()
		{
			var dense = Matrix.CreateFromArray(SparseRectangular10by5.Matrix);
			var computed = CsrMatrix.CreateFromDense(dense);
			var expected = CsrMatrix.CreateFromArrays(SparseRectangular10by5.NumRows, SparseRectangular10by5.NumCols,
				SparseRectangular10by5.CsrValues, SparseRectangular10by5.CsrColIndices, SparseRectangular10by5.CsrRowOffsets,
				true);
			comparer.AssertEqual(expected, computed);
		}

		[Fact]
        private static void TestEquality()
        {
            var full = Matrix.CreateFromArray(SparseRectangular10by5.Matrix);
            var csr = CsrMatrix.CreateFromArrays(SparseRectangular10by5.NumRows, SparseRectangular10by5.NumCols,
                SparseRectangular10by5.CsrValues, SparseRectangular10by5.CsrColIndices, SparseRectangular10by5.CsrRowOffsets,
                true);
            Assert.True(csr.Equals(full));
        }

		[SkippableFact(typeof(PerformanceBottleneckException))]
		private static void TestGetColumn()
        {
			var matrix = CsrMatrix.CreateFromArrays(SparseRectangular10by5.NumRows, SparseRectangular10by5.NumCols,
                SparseRectangular10by5.CsrValues, SparseRectangular10by5.CsrColIndices, SparseRectangular10by5.CsrRowOffsets,
                true);
            for (int j = 0; j < SparseRectangular10by5.NumCols; ++j)
            {
                Vector colExpected = DenseStrategies.GetColumn(matrix, j);
                Vector colComputed = matrix.GetColumn(j);
                comparer.AssertEqual(colExpected, colComputed);
            }
        }

        [SkippableFact(typeof(PerformanceBottleneckException))]
        private static void TestGetRow()
        {
            var matrix = CsrMatrix.CreateFromArrays(SparseRectangular10by5.NumRows, SparseRectangular10by5.NumCols,
                SparseRectangular10by5.CsrValues, SparseRectangular10by5.CsrColIndices, SparseRectangular10by5.CsrRowOffsets,
                true);
            for (int i = 0; i < SparseRectangular10by5.NumRows; ++i)
            {
                Vector rowExpected = DenseStrategies.GetRow(matrix, i);
                Vector rowComputed = matrix.GetRow(i);
                comparer.AssertEqual(rowExpected, rowComputed);
            }
        }

		[Fact]
		private static void TestGetDiagonal()
		{
			var denseMatrix = Matrix.CreateFromArray(SquareInvertible10by10.Matrix);
			var csrMatrix = CsrMatrix.CreateFromDense(denseMatrix);
			var diagonalExpected = Vector.CreateFromArray(SquareInvertible10by10.Diagonal, true);

			var diagonalComputed = Vector.CreateFromArray(csrMatrix.GetDiagonalAsArray());
			comparer.AssertEqual(diagonalExpected, diagonalComputed);

		}

		[Fact]
        private static void TestMatrixCopy()
        {
            var full = Matrix.CreateFromArray(SparseRectangular10by5.Matrix);
            var csr = CsrMatrix.CreateFromArrays(SparseRectangular10by5.NumRows, SparseRectangular10by5.NumCols,
                SparseRectangular10by5.CsrValues, SparseRectangular10by5.CsrColIndices, SparseRectangular10by5.CsrRowOffsets, 
                true);
            comparer.AssertEqual(full, csr.CopyToFullMatrix());
        }

        [Theory]
        [MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
        private static void TestMatrixMatrixMultiplication(LinearAlgebraProviderChoice providers)
        {
            TestSettings.RunMultiproviderTest(providers, delegate ()
            {
                var matrix5x5 = Matrix.CreateFromArray(SquareInvertible10by10.Matrix).GetSubmatrix(0, 5, 0, 5); //TODO: add a 5x5 denseMatrix and its products
                var matrix10x10 = Matrix.CreateFromArray(SquareInvertible10by10.Matrix);
                var ATimesMatrix5x5 = Matrix.CreateFromArray(
                    MatrixOperations.MatrixTimesMatrix(SparseRectangular10by5.Matrix, matrix5x5.CopyToArray2D()));
                var ATimesTransposeMatrix5x5 = Matrix.CreateFromArray(
                    MatrixOperations.MatrixTimesMatrix(SparseRectangular10by5.Matrix, matrix5x5.Transpose().CopyToArray2D()));
                var transposeATimesMatrix10x10 = Matrix.CreateFromArray(MatrixOperations.MatrixTimesMatrix(
                    MatrixOperations.Transpose(SparseRectangular10by5.Matrix), matrix10x10.CopyToArray2D()));
                var transposeATimesTransposeMatrix10x10 = Matrix.CreateFromArray(MatrixOperations.MatrixTimesMatrix(
                    MatrixOperations.Transpose(SparseRectangular10by5.Matrix), matrix10x10.Transpose().CopyToArray2D()));
                var matrix10x10TimesA = Matrix.CreateFromArray(
                    MatrixOperations.MatrixTimesMatrix(matrix10x10.CopyToArray2D(), SparseRectangular10by5.Matrix));
                var transposeMatrix10x10TimesA = Matrix.CreateFromArray(
                    MatrixOperations.MatrixTimesMatrix(matrix10x10.Transpose().CopyToArray2D(), SparseRectangular10by5.Matrix));
                var matrix5x5TimesTransposeA = Matrix.CreateFromArray(
                    MatrixOperations.MatrixTimesMatrix(matrix5x5.CopyToArray2D(),
                    MatrixOperations.Transpose(SparseRectangular10by5.Matrix)));
                var transposeMatrix5x5TimesTransposeA = Matrix.CreateFromArray(
                    MatrixOperations.MatrixTimesMatrix(matrix5x5.Transpose().CopyToArray2D(),
                    MatrixOperations.Transpose(SparseRectangular10by5.Matrix)));

                var A = CsrMatrix.CreateFromArrays(SparseRectangular10by5.NumRows, SparseRectangular10by5.NumCols,
                   SparseRectangular10by5.CsrValues, SparseRectangular10by5.CsrColIndices, SparseRectangular10by5.CsrRowOffsets,
                   true);

                // MultiplyRight()
                comparer.AssertEqual(ATimesMatrix5x5, A.MultiplyRight(matrix5x5, false, false));
                comparer.AssertEqual(ATimesTransposeMatrix5x5, A.MultiplyRight(matrix5x5, false, true));
                comparer.AssertEqual(transposeATimesMatrix10x10, A.MultiplyRight(matrix10x10, true, false));
                comparer.AssertEqual(transposeATimesTransposeMatrix10x10, A.MultiplyRight(matrix10x10, true, true));

                // MultiplyLeft()
                comparer.AssertEqual(matrix10x10TimesA, A.MultiplyLeft(matrix10x10, false, false));
                comparer.AssertEqual(transposeMatrix10x10TimesA, A.MultiplyLeft(matrix10x10, false, true));
                comparer.AssertEqual(matrix5x5TimesTransposeA, A.MultiplyLeft(matrix5x5, true, false));
                comparer.AssertEqual(transposeMatrix5x5TimesTransposeA, A.MultiplyLeft(matrix5x5, true, true));
            });
        }

        [Theory]
        [MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
        private static void TestMatrixVectorMultiplication(LinearAlgebraProviderChoice providers)
        {
            TestSettings.RunMultiproviderTest(providers, delegate ()
            {
                // MultiplyRight() - untransposed - rows > cols
                var A = CsrMatrix.CreateFromArrays(SparseRectangular10by5.NumRows, SparseRectangular10by5.NumCols,
                    SparseRectangular10by5.CsrValues, SparseRectangular10by5.CsrColIndices, SparseRectangular10by5.CsrRowOffsets,
                    true);
                var x5 = Vector.CreateFromArray(SparseRectangular10by5.Lhs5);
                var b10Expected = Vector.CreateFromArray(SparseRectangular10by5.Rhs10);
                Vector b10Computed = A.Multiply(x5, false);
                comparer.AssertEqual(b10Expected, b10Computed);

                // MultiplyRight() - transposed - rows > cols
                var x10 = Vector.CreateFromArray(SparseRectangular10by5.Lhs10);
                var b5Expected = Vector.CreateFromArray(SparseRectangular10by5.Rhs5);
                Vector b5Computed = A.Multiply(x10, true);
                comparer.AssertEqual(b5Expected, b5Computed);

                // MultiplyVectorSection()
                var x15 = Vector.CreateWithValue(5, 1000.0).Append(x5);
                var b20Expected = Vector.CreateWithValue(10, 10.0).Append(A.Multiply(x5, false));
                var b20Computed = Vector.CreateWithValue(20, 10.0);
                A.MultiplyVectorSection(x15, 5, b20Computed, 10);
                comparer.AssertEqual(b20Expected, b20Computed);

                // MultiplyRight() - untransposed - rows < cols
                var B = CsrMatrix.CreateFromArrays(SparseRectangular10by5.NumCols, SparseRectangular10by5.NumRows,
                    SparseRectangular10by5.CscValues, SparseRectangular10by5.CscRowIndices, SparseRectangular10by5.CscColOffsets,
                    true);
                var y10 = Vector.CreateFromArray(SparseRectangular10by5.Lhs10);
                var c5Expected = Vector.CreateFromArray(SparseRectangular10by5.Rhs5);
                Vector c5Computed = B.Multiply(y10, false);
                comparer.AssertEqual(c5Expected, c5Computed);

                // MultiplyRight() - transposed - rows < cols
                var y5 = Vector.CreateFromArray(SparseRectangular10by5.Lhs5);
                var c10Expected = Vector.CreateFromArray(SparseRectangular10by5.Rhs10);
                Vector c10Computed = B.Multiply(y5, true);
                comparer.AssertEqual(c10Expected, c10Computed);
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

                // MultiplyIntoResult() - untransposed 
                var A = CsrMatrix.CreateFromArrays(SparseRectangular10by5.NumRows, SparseRectangular10by5.NumCols,
                    SparseRectangular10by5.CsrValues, SparseRectangular10by5.CsrColIndices, SparseRectangular10by5.CsrRowOffsets,
                    true);
                var x5 = Vector.CreateFromArray(SparseRectangular10by5.Lhs5);
                var b10Expected = Vector.CreateFromArray(SparseRectangular10by5.Rhs10);
                Vector b10Computed = Vector.CreateWithValue(SparseRectangular10by5.NumRows, 1.0);
                //Vector bComputed = Vector.CreateZero(SparseRectangular10by5.NumRows);
                A.MultiplyIntoResult(x5, b10Computed, false);
                comparer.AssertEqual(b10Expected, b10Computed);

                // MultiplyIntoResult() - transposed
                var x10 = Vector.CreateFromArray(SparseRectangular10by5.Lhs10);
                var b5Expected = Vector.CreateFromArray(SparseRectangular10by5.Rhs5);
                Vector b5Computed = Vector.CreateWithValue(SparseRectangular10by5.NumCols, 1.0);
                A.MultiplyIntoResult(x10, b5Computed, true);
                comparer.AssertEqual(b5Expected, b5Computed);
            });
        }

        [Fact] //TODO: If the explicit transposition becomes abstracted behind a provider, then this should also be a Theory
        private static void TestTransposition()
        {
            var matrix = CsrMatrix.CreateFromArrays(SparseRectangular10by5.NumRows, SparseRectangular10by5.NumCols,
                SparseRectangular10by5.CsrValues, SparseRectangular10by5.CsrColIndices, SparseRectangular10by5.CsrRowOffsets,
                true);
            var transposeExpected = Matrix.CreateFromArray(MatrixOperations.Transpose(SparseRectangular10by5.Matrix));

            // TransposeToCSC()
            CscMatrix transposeCsc = matrix.TransposeToCSC(true);
            comparer.AssertEqual(transposeExpected, transposeCsc);

            // TransposeToCSR()
            CsrMatrix transposeCsr = matrix.TransposeToCSR();
            comparer.AssertEqual(transposeExpected, transposeCsr);
        }
    }
}
