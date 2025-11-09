using System;
using System.IO;
using System.Runtime.Serialization.Formatters.Binary;
using MGroup.LinearAlgebra.Commons;
using MGroup.LinearAlgebra.Exceptions;
using MGroup.LinearAlgebra.Matrices;
using MGroup.LinearAlgebra.Tests.Mocking;
using MGroup.LinearAlgebra.Tests.TestData;
using MGroup.LinearAlgebra.Tests.Utilities;
using MGroup.LinearAlgebra.Vectors;
using Xunit;

namespace MGroup.LinearAlgebra.Tests.Matrices
{
	/// <summary>
	/// Tests for <see cref="DefaultMatrix"/>.
	/// </summary>
	public static class DefaultMatrixTests
	{
		private static readonly MatrixComparer comparer = new MatrixComparer(1E-13);

		[Fact]
		internal static void TestAddition()
		{
			var A1 = new MockDefaultMatrix(SquareSingular10by10.Matrix);
			var A2 = new MockDefaultMatrix(SymmPosDef10by10.Matrix);
			var expected = Matrix.CreateFromArray(
				MatrixOperations.LinearCombination(1.0, SquareSingular10by10.Matrix, 1.0, SymmPosDef10by10.Matrix));

			// Add()
			comparer.AssertEqual(expected, A1.Add(A2));

			// AddIntoThis()
			var computed = new MockDefaultMatrix(SquareSingular10by10.Matrix);
			computed.AddIntoThis(A2);
			comparer.AssertEqual(expected, computed);
		}

		[Fact]
		internal static void TestAxpy()
		{
			var A1 = new MockDefaultMatrix(SquareSingular10by10.Matrix);
			var A2 = new MockDefaultMatrix(SymmPosDef10by10.Matrix);
			var expected = Matrix.CreateFromArray(
				MatrixOperations.LinearCombination(1.0, SquareSingular10by10.Matrix, -5.0, SymmPosDef10by10.Matrix));

			// Axpy()
			comparer.AssertEqual(expected, A1.Axpy(A2, -5.0));

			// AxpyIntoThis()
			var computed = new MockDefaultMatrix(SquareSingular10by10.Matrix);
			computed.AxpyIntoThis(A2, -5.0);
			comparer.AssertEqual(expected, computed);
		}

		[Fact]
		internal static void TestClear()
		{
			var zero = Matrix.CreateZero(SparseRectangular10by5.NumRows, SparseRectangular10by5.NumCols);
			var matrix = new MockDefaultMatrix(SparseRectangular10by5.Matrix);
			matrix.Clear();
			comparer.AssertEqual(zero, matrix);
		}

		[Fact]
		internal static void TestCopy()
		{
			var original = new MockDefaultMatrix(SparsePosDef10by10.Matrix);
			var expected = Matrix.CreateFromArray(SparsePosDef10by10.Matrix);

			// Copy(false)
			IMatrix clone1 = original.Copy(false);
			comparer.AssertEqual(expected, clone1);

			// Copy(true)
			IMatrix clone2 = original.Copy(true);
			comparer.AssertEqual(expected, clone2);

			// CopyToFullMatrix()
			Matrix clone3 = original.CopyToFullMatrix();
			comparer.AssertEqual(expected, clone3);
		}

		[Fact]
		internal static void TestDoEntrywise()
		{
			var A1 = new MockDefaultMatrix(SquareSingular10by10.Matrix);
			var A2 = new MockDefaultMatrix(SymmPosDef10by10.Matrix);
			var expected = Matrix.CreateFromArray(
				MatrixOperations.LinearCombination(3.3, SquareSingular10by10.Matrix, -5.0, SymmPosDef10by10.Matrix));

			// DoEntrywise()
			comparer.AssertEqual(expected, A1.DoEntrywise(A2, (x, y) => 3.3 * x - 5.0 * y));

			// DoEntrywiseIntoThis()
			var computed = new MockDefaultMatrix(SquareSingular10by10.Matrix);
			computed.DoEntrywiseIntoThis(A2, (x, y) => 3.3 * x - 5.0 * y);
			comparer.AssertEqual(expected, computed);
		}

		[Fact]
		internal static void TestDoToAllEntries()
		{
			var matrix = new MockDefaultMatrix(SquareSingular10by10.Matrix);

			var expected = Matrix.CreateFromArray(MatrixOperations.Round(
					MatrixOperations.Scale(Math.PI, SquareSingular10by10.Matrix),
					3)
			);

			// DoToAllEntries()
			comparer.AssertEqual(expected, matrix.DoToAllEntries(x => Math.Round(x * Math.PI, 3)));

			// DoEntrywiseIntoThis()
			var computed = new MockDefaultMatrix(SquareSingular10by10.Matrix);
			computed.DoToAllEntriesIntoThis(x => Math.Round(x * Math.PI, 3));
			comparer.AssertEqual(expected, computed);
		}

		[Fact]
		internal static void TestEquals()
		{
			// Equals(SkylineMatrix)
			var default1 = new MockDefaultMatrix(SparsePosDef10by10.Matrix);
			var skyline1 = SkylineMatrix.CreateFromArrays(SparsePosDef10by10.Order,
				SparsePosDef10by10.SkylineValues, SparsePosDef10by10.SkylineDiagOffsets, true, true);
			Assert.True(default1.Equals(skyline1));
			Assert.True(skyline1.Equals(default1));

			// Equals(CsrMatrix)
			var default2 = new MockDefaultMatrix(SparseRectangular10by5.Matrix);
			var csr2 = CsrMatrix.CreateFromArrays(SparseRectangular10by5.NumRows, SparseRectangular10by5.NumCols,
				SparseRectangular10by5.CsrValues, SparseRectangular10by5.CsrColIndices, SparseRectangular10by5.CsrRowOffsets,
				true);
			Assert.True(default2.Equals(csr2));
			Assert.True(csr2.Equals(default2));

			// Equals(CscMatrix)
			var default3 = new MockDefaultMatrix(SparseRectangular10by5.Matrix);
			var csc3 = CscMatrix.CreateFromArrays(SparseRectangular10by5.NumRows, SparseRectangular10by5.NumCols,
				SparseRectangular10by5.CscValues, SparseRectangular10by5.CscRowIndices, SparseRectangular10by5.CscColOffsets,
				true);
			Assert.True(default3.Equals(csc3));
			Assert.True(csc3.Equals(default3));

			// Equals(Matrix)
			var default4 = new MockDefaultMatrix(SparsePosDef10by10.Matrix);
			var full4 = Matrix.CreateFromArray(SparsePosDef10by10.Matrix);
			Assert.True(default4.Equals(full4));
			Assert.True(full4.Equals(default4));
		}

		[Fact]
		internal static void TestGetDiagonal()
		{
			var matrix = new MockDefaultMatrix(SparsePosDef10by10.Matrix);

			// GetDiagonal()
			comparer.AssertEqual(Vector.CreateFromArray(SparsePosDef10by10.Diagonal), matrix.GetDiagonal());

			// GetDiagonalAsArray()
			comparer.AssertEqual(SparsePosDef10by10.Diagonal, matrix.GetDiagonalAsArray());
		}


		[SkippableFact(typeof(PerformanceBottleneckException))]
		internal static void TestGetColumn()
		{
			var matrix = new MockDefaultMatrix(RectangularFullRank10by5.Matrix);
			for (int j = 0; j < RectangularFullRank10by5.NumCols; ++j)
			{
				Vector colExpected = DenseStrategies.GetColumn(matrix, j);
				Vector colComputed = matrix.GetColumn(j);
				comparer.AssertEqual(colExpected, colComputed);
			}
		}

		[SkippableFact(typeof(PerformanceBottleneckException))]
		internal static void TestGetRow()
		{
			var matrix = new MockDefaultMatrix(RectangularFullRank10by5.Matrix);
			for (int i = 0; i < RectangularFullRank10by5.NumRows; ++i)
			{
				Vector rowExpected = DenseStrategies.GetRow(matrix, i);
				Vector rowComputed = matrix.GetRow(i);
				comparer.AssertEqual(rowExpected, rowComputed);
			}
		}

		[Fact]
		internal static void TestGetSubmatrix()
		{
			var matrix = new MockDefaultMatrix(SparsePosDef10by10.Matrix);
			var expected = Matrix.CreateFromArray(SparsePosDef10by10.SubmatrixRows13Cols45);

			// GetSubmatrix(index arrays)
			IMatrix submatrixA = matrix.GetSubmatrix(new int[] { 1, 2, 3 }, new int[] { 4, 5 });
			comparer.AssertEqual(expected, submatrixA);

			// GetSubmatrix(from/to indices)
			IMatrix submatrixB = matrix.GetSubmatrix(1, 4, 4, 6);
			comparer.AssertEqual(expected, submatrixB);
		}

		[Fact]
		internal static void TestLinearCombination()
		{
			var A1 = new MockDefaultMatrix(SquareSingular10by10.Matrix);
			double scalar1 = 2.0;
			var A2 = new MockDefaultMatrix(SymmPosDef10by10.Matrix);
			double scalar2 = 3.5;
			var expected = Matrix.CreateFromArray(
				MatrixOperations.LinearCombination(scalar1, SquareSingular10by10.Matrix, scalar2, SymmPosDef10by10.Matrix));

			// LinearCombination()
			comparer.AssertEqual(expected, A1.LinearCombination(scalar1, A2, scalar2));

			// LinearCombinationIntoThis()
			var computed = new MockDefaultMatrix(A1.CopyToArray2D());
			computed.LinearCombinationIntoThis(scalar1, A2, scalar2);
			comparer.AssertEqual(expected, computed);
		}

		[Fact]
		internal static void TestMatrixMatrixMultiplication()
		{
			var A1 = new MockDefaultMatrix(SquareSingular10by10.Matrix);
			var A2 = new MockDefaultMatrix(RectangularFullRank10by5.Matrix);
			var expectedA1TimesA2 = Matrix.CreateFromArray(
				MatrixOperations.MatrixTimesMatrix(SquareSingular10by10.Matrix, RectangularFullRank10by5.Matrix));
			var expectedTransposeA2TimesA1 = Matrix.CreateFromArray(
				MatrixOperations.MatrixTimesMatrix(
					MatrixOperations.Transpose(RectangularFullRank10by5.Matrix), SquareSingular10by10.Matrix));
			var expectedTransposeA2TimesTransposeA1 = Matrix.CreateFromArray(
				MatrixOperations.MatrixTimesMatrix(
					MatrixOperations.Transpose(RectangularFullRank10by5.Matrix),
					MatrixOperations.Transpose(SquareSingular10by10.Matrix)
				)
			);

			// MultiplyRight() without transposition
			comparer.AssertEqual(expectedA1TimesA2, A1.MultiplyRight(A2, false, false));

			// MultiplyRight() with 1 transposition
			comparer.AssertEqual(expectedTransposeA2TimesA1, A2.MultiplyRight(A1, true, false));

			// MultiplyRight() with 2 transpositions
			comparer.AssertEqual(expectedTransposeA2TimesTransposeA1, A2.MultiplyRight(A1, true, true));

			// MultiplyRight() with incorrect dimensions
			Assert.Throws<NonMatchingDimensionsException>(() => A2.MultiplyRight(A1, false, false));

			// MultiplyLeft without transposition
			comparer.AssertEqual(expectedA1TimesA2, A2.MultiplyLeft(A1));

			// MultiplyLeft with 1 transposition
			comparer.AssertEqual(expectedTransposeA2TimesA1, A1.MultiplyLeft(A2, false, true));

			// MultiplyLeft() with 2 transpositions
			comparer.AssertEqual(expectedTransposeA2TimesTransposeA1, A1.MultiplyLeft(A2, true, true));

			// MultiplyLeft() with incorrect dimensions
			Assert.Throws<NonMatchingDimensionsException>(() => A1.MultiplyLeft(A2, false, false));
		}

		[Fact]
		internal static void TestMatrixVectorMultiplication()
		{
			var mvChecker = new MatrixDenseVectorMultiplicationChecker(
				(A, x, tranpose) => (Vector)(A.Multiply(x, tranpose)),
				(A, x, y, tranpose) => A.MultiplyIntoResult(x, y, tranpose)
			);
			mvChecker.Tolerance = 1E-13;

			// rectangular 10-by-5
			var A1 = new MockDefaultMatrix(RectangularFullRank10by5.Matrix);
			mvChecker.CheckAllMultiplications(A1, RectangularFullRank10by5.Lhs5, RectangularFullRank10by5.Rhs10, false);

			// rectangular 5-by-10
			mvChecker.CheckAllMultiplications(A1, RectangularFullRank10by5.Lhs10, RectangularFullRank10by5.Rhs5, true);

			// square invertible 10-by-10
			var A3 = new MockDefaultMatrix(SquareInvertible10by10.Matrix);
			mvChecker.CheckAllMultiplications(A3, SquareInvertible10by10.Lhs, SquareInvertible10by10.Rhs, false);

			// square singular 10-by-10 (rank = 8)
			var A4 = new MockDefaultMatrix(SquareSingular10by10.Matrix);
			mvChecker.CheckAllMultiplications(A4, SquareSingular10by10.Lhs, SquareSingular10by10.Rhs, false);

			// square singular 10-by-10 (rank = 9)
			var A5 = new MockDefaultMatrix(SquareSingularSingleDeficiency10by10.Matrix);
			mvChecker.CheckAllMultiplications(A5, SquareSingularSingleDeficiency10by10.Lhs,
				SquareSingularSingleDeficiency10by10.Rhs, false);
		}

		[Fact]
		internal static void TestMatrixVectorMultiplicationIntoResult()
		{
			// The result vectors will first be set to some non zero values to make sure that the result overwrites 
			// them instead of being added to them.

			// MultiplyIntoResult() - untransposed 
			var A1 = new MockDefaultMatrix(RectangularFullRank10by5.Matrix);
			var x1 = Vector.CreateFromArray(RectangularFullRank10by5.Lhs5);
			var b1Expected = Vector.CreateFromArray(RectangularFullRank10by5.Rhs10);
			Vector b1Computed = Vector.CreateWithValue(A1.NumRows, 1.0);
			A1.MultiplyIntoResult(x1, b1Computed, false);
			comparer.AssertEqual(b1Expected, b1Computed);

			// MultiplyIntoResult() - transposed
			var x2 = Vector.CreateFromArray(RectangularFullRank10by5.Lhs10);
			var b2Expected = Vector.CreateFromArray(RectangularFullRank10by5.Rhs5);
			Vector b2Computed = Vector.CreateWithValue(A1.NumColumns, 1.0);
			A1.MultiplyIntoResult(x2, b2Computed, true);
			comparer.AssertEqual(b2Expected, b2Computed);
		}

		[Fact]
		internal static void TestReduce()
		{
			var matrix = new MockDefaultMatrix(RectangularFullRank10by5.Matrix);
			double maxAbsExpected = MatrixOperations.ReduceMaxAbs(RectangularFullRank10by5.Matrix);

			double maxAbsComputed = matrix.Reduce(0.0,
				(x, amax) => Math.Abs(x) > amax ? Math.Abs(x) : amax,
				(nz, amax) => amax,
				amax => amax);

			Assert.Equal(maxAbsExpected, maxAbsComputed);
		}

		[Fact]
		internal static void TestScale()
		{
			var matrix = new MockDefaultMatrix(RectangularFullRank10by5.Matrix);
			double scalar = 5.0;
			var expected = Matrix.CreateFromArray(MatrixOperations.Scale(scalar, RectangularFullRank10by5.Matrix));

			// Scale()
			comparer.AssertEqual(expected, matrix.Scale(scalar));

			// ScaleIntoThis()
			var computed = new MockDefaultMatrix(matrix.CopyToArray2D());
			computed.ScaleIntoThis(scalar);
			comparer.AssertEqual(expected, computed);
		}

		[Fact]
		internal static void TestSerialization()
		{
			var originalMatrix = new MockDefaultMatrix(RectangularFullRank10by5.Matrix);
			var formatter = new BinaryFormatter();
			using (var stream = new MemoryStream())
			{
				formatter.Serialize(stream, originalMatrix);
				stream.Seek(0, SeekOrigin.Begin);
				var deserializedMatrix = (MockDefaultMatrix)formatter.Deserialize(stream);

				Assert.True(originalMatrix.Equals(deserializedMatrix));
			}
		}

		[Fact]
		internal static void TestSubtraction()
		{
			var A1 = new MockDefaultMatrix(SquareSingular10by10.Matrix);
			var A2 = new MockDefaultMatrix(SymmPosDef10by10.Matrix);
			var expected = Matrix.CreateFromArray(
				MatrixOperations.LinearCombination(1.0, SquareSingular10by10.Matrix, -1.0, SymmPosDef10by10.Matrix));

			// Subtract()
			comparer.AssertEqual(expected, A1.Subtract(A2));

			// SubtractIntoThis()
			var computed = new MockDefaultMatrix(A1.CopyToArray2D());
			computed.SubtractIntoThis(A2);
			comparer.AssertEqual(expected, computed);
		}

		[Fact]
		internal static void TestSet()
		{
			var computed = new MockDefaultMatrix(SquareSingular10by10.Matrix);
			computed.Set(0, 0, -1000);
			computed.Set(6, 2, 500.39);
			computed.Set(1, 8, 0.2222333444);

			double[,] expected = SquareSingular10by10.Matrix;
			expected[0, 0] = -1000;
			expected[6, 2] = 500.39;
			expected[1, 8] = 0.2222333444;

			comparer.AssertEqual(Matrix.CreateFromArray(expected), computed);
		}

		[Fact]
		internal static void TestTransposition()
		{
			// Square
			var A1 = new MockDefaultMatrix(SquareSingular10by10.Matrix);
			var A1TransposeExpected = MatrixOperations.Transpose(SquareSingular10by10.Matrix);
			IMatrix A1TransposeComputed = A1.Transpose();
			comparer.AssertEqual(A1TransposeExpected, A1TransposeComputed.CopyToArray2D());

			// Rectangular
			var A2 = new MockDefaultMatrix(RectangularFullRank10by5.Matrix);
			var A2TransposeExpected = MatrixOperations.Transpose(RectangularFullRank10by5.Matrix);
			IMatrix A2TransposeComputed = A2.Transpose();
			comparer.AssertEqual(A2TransposeExpected, A2TransposeComputed.CopyToArray2D());
		}
	}
}
