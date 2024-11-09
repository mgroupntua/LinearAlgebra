using System;
using System.IO;
using System.Linq;
using System.Runtime.Serialization.Formatters.Binary;
using MGroup.LinearAlgebra.Commons;
using MGroup.LinearAlgebra.Exceptions;
using MGroup.LinearAlgebra.Implementations;
using MGroup.LinearAlgebra.Matrices;
using MGroup.LinearAlgebra.Tests.TestData;
using MGroup.LinearAlgebra.Tests.Utilities;
using MGroup.LinearAlgebra.Vectors;
using Xunit;

namespace MGroup.LinearAlgebra.Tests.Matrices
{
	/// <summary>
	/// Tests for <see cref="Matrix"/>.
	/// Authors: Serafeim Bakalakos
	/// </summary>
	public static class MatrixTests
	{
		private static readonly MatrixComparer comparer = new MatrixComparer(1E-13);

		[Theory]
		[MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
		internal static void TestAddition(IImplementationProvider provider)
		{
			TestSettings.RunMultiproviderTest(provider, delegate ()
			{
				var A1 = Matrix.CreateFromArray(SquareSingular10by10.Matrix);
				var A2 = Matrix.CreateFromArray(SymmPosDef10by10.Matrix);
				var expected = Matrix.CreateFromArray(
					MatrixOperations.LinearCombination(1.0, SquareSingular10by10.Matrix, 1.0, SymmPosDef10by10.Matrix));

				// operator+
				comparer.AssertEqual(expected, A1 + A2);
			});
		}

		[Theory]
		[MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
		internal static void TestAxpyColumn(IImplementationProvider provider)
		{
			TestSettings.RunMultiproviderTest(provider, delegate ()
			{
				var A1 = Matrix.CreateFromArray(SquareSingular10by10.Matrix);
				var columnVector = Vector.CreateWithValue(10, 1.0);
				int colIdx = 3;
				double coeff = 2.0;

				double[,] expected = SquareSingular10by10.Matrix;
				AxpyColumn(expected, colIdx, columnVector.RawData, coeff);
				A1.AxpyColumn(colIdx, coeff, columnVector);

				comparer.AssertEqual(Matrix.CreateFromArray(expected), A1);
			});
		}

		[Fact]
		private static void TestClear()
		{
			var zero = Matrix.CreateZero(SparseRectangular10by5.NumRows, SparseRectangular10by5.NumCols);
			var matrix = Matrix.CreateFromArray(SparseRectangular10by5.Matrix);
			matrix.Clear();
			comparer.AssertEqual(zero, matrix);
		}

		[Fact]
		private static void TestCreateWithValue()
		{
			int m = 35;
			int n = 20;
			double val = -3.33;
			var matrix = Matrix.CreateWithValue(m, n, val);

			var expected = new double[m, n];
			LinearAlgebra.Tests.Utilities.ArrayUtilities.SetAll(expected, val);
			var comparer = new MatrixComparer(1E-15);
			comparer.AssertEqual(expected, matrix);
		}

		[Fact]
		private static void TestEquality()
		{
			// Equals(SkylineMatrix)
			var full1 = Matrix.CreateFromArray(SparsePosDef10by10.Matrix);
			var skyline1 = SkylineMatrix.CreateFromArrays(SparsePosDef10by10.Order,
				SparsePosDef10by10.SkylineValues, SparsePosDef10by10.SkylineDiagOffsets, true, true);
			Assert.True(full1.Equals(skyline1));

			// Equals(CsrMatrix)
			var full2 = Matrix.CreateFromArray(SparseRectangular10by5.Matrix);
			var csr2 = CsrMatrix.CreateFromArrays(SparseRectangular10by5.NumRows, SparseRectangular10by5.NumCols,
				SparseRectangular10by5.CsrValues, SparseRectangular10by5.CsrColIndices, SparseRectangular10by5.CsrRowOffsets,
				true);
			Assert.True(full2.Equals(csr2));

			// Equals(CscMatrix)
			var full3 = Matrix.CreateFromArray(SparseRectangular10by5.Matrix);
			var csc3 = CscMatrix.CreateFromArrays(SparseRectangular10by5.NumRows, SparseRectangular10by5.NumCols,
				SparseRectangular10by5.CscValues, SparseRectangular10by5.CscRowIndices, SparseRectangular10by5.CscColOffsets,
				true);
			Assert.True(full3.Equals(csc3));
		}

		[SkippableFact(typeof(PerformanceBottleneckException))]
		private static void TestGetColumn()
		{
			var matrix = Matrix.CreateFromArray(RectangularFullRank10by5.Matrix);
			for (int j = 0; j < RectangularFullRank10by5.NumCols; ++j)
			{
				Vector colExpected = DenseStrategies.GetColumn(matrix, j);
				Vector colComputed = matrix.GetColumn(j);
				comparer.AssertEqual(colExpected, colComputed);
			}
		}

		[SkippableFact(typeof(PerformanceBottleneckException))]
		private static void TestGetRow()
		{
			var matrix = Matrix.CreateFromArray(RectangularFullRank10by5.Matrix);
			for (int i = 0; i < RectangularFullRank10by5.NumRows; ++i)
			{
				Vector rowExpected = DenseStrategies.GetRow(matrix, i);
				Vector rowComputed = matrix.GetRow(i);
				comparer.AssertEqual(rowExpected, rowComputed);
			}
		}

		[Theory]
		[MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
		internal static void TestInvertAndDeterminant(IImplementationProvider provider)
		{
			TestSettings.RunMultiproviderTest(provider, delegate ()
			{
				var matrix = Matrix.CreateFromArray(new double[,]
				{
					{ 87.5, 0, 0, 1 },
					{ 90, 0, 0, 1 },
					{ 90, 2.5, 225, 1 },
					{ 87.5, 2.5, 218.75, 1 }
				});
				(Matrix inverse, double determinant) = matrix.InvertAndDeterminant();

				double detExpected = 39.0625;
				var inverseExpected = Matrix.CreateFromArray(new double[,]
				{
					{ -0.4, 0.4, 0, 0 },
					{ -14.4, 14, -14, 14.4 },
					{ 0.16, -0.16, 0.16, -0.16 },
					{ 36, -35, 0, 0 }
				});

				// operator+
				Assert.Equal(detExpected, determinant, 9);
				comparer.AssertEqual(inverseExpected, inverse);
			});
		}

		[Theory]
		[MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
		private static void TestLinearCombination(IImplementationProvider provider)
		{
			TestSettings.RunMultiproviderTest(provider, delegate ()
			{
				var A1 = Matrix.CreateFromArray(SquareSingular10by10.Matrix);
				double scalar1 = 2.0;
				var A2 = Matrix.CreateFromArray(SymmPosDef10by10.Matrix);
				double scalar2 = 3.5;
				var expected = Matrix.CreateFromArray(
					MatrixOperations.LinearCombination(scalar1, SquareSingular10by10.Matrix, scalar2, SymmPosDef10by10.Matrix));

				// LinearCombination()
				comparer.AssertEqual(expected, A1.LinearCombination(scalar1, A2, scalar2));

				// LinearCombinationIntoThis()
				Matrix temp = A1.Copy();
				temp.LinearCombinationIntoThis(scalar1, A2, scalar2);
				comparer.AssertEqual(expected, temp);
			});
		}

		[Theory]
		[MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
		private static void TestMatrixMatrixMultiplication(IImplementationProvider provider)
		{
			TestSettings.RunMultiproviderTest(provider, delegate ()
			{
				var A1 = Matrix.CreateFromArray(SquareSingular10by10.Matrix);
				var A2 = Matrix.CreateFromArray(RectangularFullRank10by5.Matrix);
				var expectedA1TimesA2 = Matrix.CreateFromArray(
					MatrixOperations.MatrixTimesMatrix(SquareSingular10by10.Matrix, RectangularFullRank10by5.Matrix));
				var expectedTransposeA2TimesA1 = Matrix.CreateFromArray(
					MatrixOperations.MatrixTimesMatrix(
						MatrixOperations.Transpose(RectangularFullRank10by5.Matrix), SquareSingular10by10.Matrix));

				// MultiplyRight() without transposition
				comparer.AssertEqual(expectedA1TimesA2, A1.MultiplyRight(A2, false, false));

				// operator*
				comparer.AssertEqual(expectedA1TimesA2, A1 * A2);

				// MultiplyRight() with transposition
				comparer.AssertEqual(expectedTransposeA2TimesA1, A2.MultiplyRight(A1, true, false));

				// MultiplyRight() with incorrect dimensions
				Assert.Throws<NonMatchingDimensionsException>(() => A2.MultiplyRight(A1, false, false));
			});
		}

		[Theory]
		[MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
		private static void TestMatrixVectorMultiplication(IImplementationProvider provider)
		{
			TestSettings.RunMultiproviderTest(provider, delegate ()
			{
				var mvChecker = new MatrixDenseVectorMultiplicationChecker(
				(A, x, tranpose) => ((Matrix)A).Multiply(x, tranpose),
				(A, x, y, tranpose) => ((Matrix)A).MultiplyIntoResult(x, y, tranpose));
				mvChecker.Tolerance = 1E-13;

				// rectangular 10-by-5
				var A1 = Matrix.CreateFromArray(RectangularFullRank10by5.Matrix);
				mvChecker.CheckAllMultiplications(A1, RectangularFullRank10by5.Lhs5, RectangularFullRank10by5.Rhs10, false);

				// rectangular 5-by-10
				mvChecker.CheckAllMultiplications(A1, RectangularFullRank10by5.Lhs10, RectangularFullRank10by5.Rhs5, true);

				// square invertible 10-by-10
				var A3 = Matrix.CreateFromArray(SquareInvertible10by10.Matrix);
				mvChecker.CheckAllMultiplications(A3, SquareInvertible10by10.Lhs, SquareInvertible10by10.Rhs, false);

				// square singular 10-by-10 (rank = 8)
				var A4 = Matrix.CreateFromArray(SquareSingular10by10.Matrix);
				mvChecker.CheckAllMultiplications(A4, SquareSingular10by10.Lhs, SquareSingular10by10.Rhs, false);

				// square singular 10-by-10 (rank = 9)
				var A5 = Matrix.CreateFromArray(SquareSingularSingleDeficiency10by10.Matrix);
				mvChecker.CheckAllMultiplications(A5, SquareSingularSingleDeficiency10by10.Lhs,
					SquareSingularSingleDeficiency10by10.Rhs, false);
			});
		}

		[Theory]
		[MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
		private static void TestMatrixVectorMultiplicationIntoResult(IImplementationProvider provider)
		{
			TestSettings.RunMultiproviderTest(provider, delegate ()
			{
				// The result vectors will first be set to some non zero values to make sure that the result overwrites 
				// them instead of being added to them.

				// MultiplyIntoResult() - untransposed 
				var A1 = Matrix.CreateFromArray(RectangularFullRank10by5.Matrix);
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
			});
		}

		[Theory]
		[MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
		private static void TestScaling(IImplementationProvider provider)
		{
			TestSettings.RunMultiproviderTest(provider, delegate ()
			{
				var matrix = Matrix.CreateFromArray(RectangularFullRank10by5.Matrix);
				double scalar = 5.0;
				var expected = Matrix.CreateFromArray(MatrixOperations.Scale(scalar, RectangularFullRank10by5.Matrix));

				// Scale()
				comparer.AssertEqual(expected, matrix.Scale(scalar));

				// ScaleIntoThis()
				Matrix temp = matrix.Copy();
				temp.ScaleIntoThis(scalar);
				comparer.AssertEqual(expected, temp);

				// operator*
				comparer.AssertEqual(expected, scalar * matrix);
			});
		}

		[Fact]
		private static void TestSerialization()
		{
			var originalMatrix = Matrix.CreateFromArray(RectangularFullRank10by5.Matrix);
			var formatter = new BinaryFormatter();
			using (var stream = new MemoryStream())
			{
				formatter.Serialize(stream, originalMatrix);
				stream.Seek(0, SeekOrigin.Begin);
				var deserializedMatrix = (Matrix)formatter.Deserialize(stream);

				Assert.True(originalMatrix.Equals(deserializedMatrix));
			}
		}

		[Fact]
		private static void TestSetAll()
		{
			int m = 35;
			int n = 20;
			double val = -3.33;
			var matrix = Matrix.CreateZero(m, n);
			matrix.SetAll(val);

			var expected = new double[m, n];
			LinearAlgebra.Tests.Utilities.ArrayUtilities.SetAll(expected, val);
			var comparer = new MatrixComparer(1E-15);
			comparer.AssertEqual(expected, matrix);
		}

		[Theory]
		[MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
		private static void TestSubtraction(IImplementationProvider provider)
		{
			TestSettings.RunMultiproviderTest(provider, delegate ()
			{
				var A1 = Matrix.CreateFromArray(SquareSingular10by10.Matrix);
				var A2 = Matrix.CreateFromArray(SymmPosDef10by10.Matrix);
				var expected = Matrix.CreateFromArray(
					MatrixOperations.LinearCombination(1.0, SquareSingular10by10.Matrix, -1.0, SymmPosDef10by10.Matrix));

				// operator+
				comparer.AssertEqual(expected, A1 - A2);
			});
		}

		[Fact]
		private static void TestTransposition()
		{
			// square
			var A1 = Matrix.CreateFromArray(SquareSingular10by10.Matrix);
			var A1TransposeExpected = MatrixOperations.Transpose(SquareSingular10by10.Matrix);
			Matrix A1TransposeComputed = A1.Transpose();
			comparer.AssertEqual(A1TransposeExpected, A1TransposeComputed.CopyToArray2D());

			// rectangular
			var A2 = Matrix.CreateFromArray(RectangularFullRank10by5.Matrix);
			var A2TransposeExpected = MatrixOperations.Transpose(RectangularFullRank10by5.Matrix);
			Matrix A2TransposeComputed = A2.Transpose();
			comparer.AssertEqual(A2TransposeExpected, A2TransposeComputed.CopyToArray2D());
		}

		private static void AxpyColumn(double[,] matrix, int colIdx, double[] columnValues, double coeff)
		{
			for (int i = 0; i < matrix.GetLength(0); i++)
			{
				matrix[i, colIdx] += coeff * columnValues[i];
			}
		}
	}
}
