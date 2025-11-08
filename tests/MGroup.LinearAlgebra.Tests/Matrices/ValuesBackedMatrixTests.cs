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
	/// Tests for <see cref="ValuesBackedMatrix{TMatrix}"/>.
	/// </summary>
	public static class ValuesBackedMatrixTests
	{
		private static readonly MatrixComparer comparer = new MatrixComparer(1E-13);

		internal static void TestAddition()
		{
			var A1 = new MockValuesBackedMatrix(
				SquareInvertible10by10.Order, SquareInvertible10by10.Order, SquareInvertible10by10.MatrixAsRowMajor);
			var A2 = new MockValuesBackedMatrix(
				SquareSingular10by10.Order, SquareSingular10by10.Order, SquareSingular10by10.MatrixAsRowMajor);
			var expected = Matrix.CreateFromArray(
				MatrixOperations.LinearCombination(1.0, SquareInvertible10by10.Matrix, 1.0, SquareSingular10by10.Matrix));

			// Add()
			comparer.AssertEqual(expected, A1.Add(A2));

			// AddIntoThis()
			var computed = new MockValuesBackedMatrix(
				SquareInvertible10by10.Order, SquareInvertible10by10.Order, SquareInvertible10by10.MatrixAsRowMajor);
			computed.AddIntoThis(A2);
			comparer.AssertEqual(expected, computed);
		}

		[Fact]
		internal static void TestAxpy()
		{
			var A1 = new MockValuesBackedMatrix(
				SquareInvertible10by10.Order, SquareInvertible10by10.Order, SquareInvertible10by10.MatrixAsRowMajor);
			var A2 = new MockValuesBackedMatrix(
				SquareSingular10by10.Order, SquareSingular10by10.Order, SquareSingular10by10.MatrixAsRowMajor);
			var expected = Matrix.CreateFromArray(
				MatrixOperations.LinearCombination(1.0, SquareInvertible10by10.Matrix, -5.0, SquareSingular10by10.Matrix));

			// Axpy()
			comparer.AssertEqual(expected, A1.Axpy(A2, -5.0));

			// AxpyIntoThis()
			var computed = new MockValuesBackedMatrix(
				SquareInvertible10by10.Order, SquareInvertible10by10.Order, SquareInvertible10by10.MatrixAsRowMajor);
			computed.AxpyIntoThis(A2, -5.0);
			comparer.AssertEqual(expected, computed);
		}

		[Fact]
		internal static void TestClear()
		{
			var zero = Matrix.CreateZero(SquareInvertible10by10.NumRows, SquareInvertible10by10.NumCols);
			var matrix = new MockValuesBackedMatrix(
				SquareInvertible10by10.Order, SquareInvertible10by10.Order, SquareInvertible10by10.MatrixAsRowMajor);
			matrix.Clear();
			comparer.AssertEqual(zero, matrix);
		}

		[Fact]
		internal static void TestCopy()
		{
			var original = new MockValuesBackedMatrix(
				SquareInvertible10by10.Order, SquareInvertible10by10.Order, SquareInvertible10by10.MatrixAsRowMajor);
			var expected = Matrix.CreateFromArray(SquareInvertible10by10.Matrix);

			// Copy(false)
			IMatrix clone1 = original.Copy(false);
			comparer.AssertEqual(expected, clone1);

			// Copy(true)
			IMatrix clone2 = original.Copy(true);
			comparer.AssertEqual(expected, clone2);
		}

		[Fact]
		internal static void TestDoEntrywise()
		{
			var A1 = new MockValuesBackedMatrix(
				SquareInvertible10by10.Order, SquareInvertible10by10.Order, SquareInvertible10by10.MatrixAsRowMajor);
			var A2 = new MockValuesBackedMatrix(
				SquareSingular10by10.Order, SquareSingular10by10.Order, SquareSingular10by10.MatrixAsRowMajor);
			var expected = Matrix.CreateFromArray(
				MatrixOperations.LinearCombination(3.3, SquareInvertible10by10.Matrix, -5.0, SquareSingular10by10.Matrix));

			// DoEntrywise()
			comparer.AssertEqual(expected, A1.DoEntrywise(A2, (x, y) => 3.3 * x - 5.0 * y));

			// DoEntrywiseIntoThis()
			var computed = new MockValuesBackedMatrix(
				SquareInvertible10by10.Order, SquareInvertible10by10.Order, SquareInvertible10by10.MatrixAsRowMajor);
			computed.DoEntrywiseIntoThis(A2, (x, y) => 3.3 * x - 5.0 * y);
			comparer.AssertEqual(expected, computed);
		}

		[Fact]
		internal static void TestDoToAllEntries()
		{
			var matrix = new MockValuesBackedMatrix(
				SquareInvertible10by10.Order, SquareInvertible10by10.Order, SquareInvertible10by10.MatrixAsRowMajor);

			var expected = Matrix.CreateFromArray(MatrixOperations.Round(
					MatrixOperations.Scale(Math.PI, SquareInvertible10by10.Matrix),
					3)
			);

			// DoToAllEntries()
			comparer.AssertEqual(expected, matrix.DoToAllEntries(x => Math.Round(x * Math.PI, 3)));

			// DoEntrywiseIntoThis()
			var computed = new MockValuesBackedMatrix(
				SquareInvertible10by10.Order, SquareInvertible10by10.Order, SquareInvertible10by10.MatrixAsRowMajor);
			computed.DoToAllEntriesIntoThis(x => Math.Round(x * Math.PI, 3));
			comparer.AssertEqual(expected, computed);
		}

		[Fact]
		internal static void TestLinearCombination()
		{
			var A1 = new MockValuesBackedMatrix(
				SquareInvertible10by10.Order, SquareInvertible10by10.Order, SquareInvertible10by10.MatrixAsRowMajor);
			var A2 = new MockValuesBackedMatrix(
				SquareSingular10by10.Order, SquareSingular10by10.Order, SquareSingular10by10.MatrixAsRowMajor);
			double scalar1 = 2.0;
			double scalar2 = 3.5;
			var expected = Matrix.CreateFromArray(
				MatrixOperations.LinearCombination(scalar1, SquareInvertible10by10.Matrix, scalar2, SquareSingular10by10.Matrix));

			// LinearCombination()
			comparer.AssertEqual(expected, A1.LinearCombination(scalar1, A2, scalar2));

			// LinearCombinationIntoThis()
			var computed = new MockValuesBackedMatrix(
				SquareInvertible10by10.Order, SquareInvertible10by10.Order, SquareInvertible10by10.MatrixAsRowMajor);
			computed.LinearCombinationIntoThis(scalar1, A2, scalar2);
			comparer.AssertEqual(expected, computed);
		}

		[Fact]
		internal static void TestScale()
		{
			var matrix = new MockValuesBackedMatrix(
				SquareInvertible10by10.Order, SquareInvertible10by10.Order, SquareInvertible10by10.MatrixAsRowMajor);
			double scalar = 5.0;
			var expected = Matrix.CreateFromArray(MatrixOperations.Scale(scalar, SquareInvertible10by10.Matrix));

			// Scale()
			comparer.AssertEqual(expected, matrix.Scale(scalar));

			// ScaleIntoThis()
			var computed = new MockValuesBackedMatrix(
				SquareInvertible10by10.Order, SquareInvertible10by10.Order, SquareInvertible10by10.MatrixAsRowMajor);
			computed.ScaleIntoThis(scalar);
			comparer.AssertEqual(expected, computed);
		}

		[Fact]
		internal static void TestSerialization()
		{
			var originalMatrix = new MockValuesBackedMatrix(
				SquareInvertible10by10.Order, SquareInvertible10by10.Order, SquareInvertible10by10.MatrixAsRowMajor);
			var formatter = new BinaryFormatter();
			using (var stream = new MemoryStream())
			{
				formatter.Serialize(stream, originalMatrix);
				stream.Seek(0, SeekOrigin.Begin);
				var deserializedMatrix = (MockValuesBackedMatrix)formatter.Deserialize(stream);

				Assert.True(originalMatrix.Equals(deserializedMatrix));
			}
		}

		[Fact]
		internal static void TestSubtraction()
		{
			var A1 = new MockValuesBackedMatrix(
				SquareInvertible10by10.Order, SquareInvertible10by10.Order, SquareInvertible10by10.MatrixAsRowMajor);
			var A2 = new MockValuesBackedMatrix(
				SquareSingular10by10.Order, SquareSingular10by10.Order, SquareSingular10by10.MatrixAsRowMajor);
			var expected = Matrix.CreateFromArray(
				MatrixOperations.LinearCombination(1.0, SquareInvertible10by10.Matrix, -1.0, SquareSingular10by10.Matrix));

			// Subtract()
			comparer.AssertEqual(expected, A1.Subtract(A2));

			// SubtractIntoThis()
			var computed = new MockValuesBackedMatrix(
				SquareInvertible10by10.Order, SquareInvertible10by10.Order, SquareInvertible10by10.MatrixAsRowMajor);
			computed.SubtractIntoThis(A2);
			comparer.AssertEqual(expected, computed);
		}
	}
}
