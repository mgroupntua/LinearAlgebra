//TODO: These should be Debug only
//TODO: when a check can be either <= (with offsets) or == (without offsets), rename the latter as ~Exact.
namespace MGroup.LinearAlgebra.Commons
{
	using System;

	using MGroup.LinearAlgebra.Exceptions;
	using MGroup.LinearAlgebra.Matrices;
	using MGroup.LinearAlgebra.Vectors;

	/// <summary>
	/// Checks for the compatibility of vector or matrix dimensions for various operations and for the indices. Also templates 
	/// for exception messages.
	/// </summary>
	public static class Preconditions
	{
		public static bool AreSameMatrixDimensions(IIndexable2D matrix1, IIndexable2D matrix2)
		{
			if ((matrix1.NumRows != matrix2.NumRows) || (matrix1.NumColumns != matrix2.NumColumns))
			{
				return false;
			}
			else
			{
				return true;
			}
		}

		public static void CheckIndex1D(IReadOnlyVector vector, int idx)
		{
			if ((idx < 0) || (idx >= vector.Length))
			{
				throw new IndexOutOfRangeException($"Cannot access index {idx} in a vector with length = {vector.Length}");
			}
		}

		public static void CheckIndexCol(IIndexable2D matrix, int colIdx)
		{
			if ((colIdx < 0) || (colIdx >= matrix.NumColumns))
			{
				throw new IndexOutOfRangeException($"Cannot access column {colIdx} in a"
					+ $" {matrix.NumRows}-by-{matrix.NumColumns} matrix");
			}
		}

		public static void CheckIndexRow(IIndexable2D matrix, int rowIdx)
		{
			if ((rowIdx < 0) || (rowIdx >= matrix.NumRows))
			{
				throw new IndexOutOfRangeException($"Cannot access row {rowIdx} in a"
					+ $" {matrix.NumRows}-by-{matrix.NumColumns} matrix");
			}
		}

		public static void CheckIndices(IIndexable2D matrix, int rowIdx, int colIdx)
		{
			if ((rowIdx < 0) || (rowIdx >= matrix.NumRows))
			{
				throw new IndexOutOfRangeException($"Cannot access row {rowIdx} in a"
					+ $" {matrix.NumRows}-by-{matrix.NumColumns} matrix");
			}

			if ((colIdx < 0) || (colIdx >= matrix.NumColumns))
			{
				throw new IndexOutOfRangeException($"Cannot access column {colIdx} in a"
					+ $" {matrix.NumRows}-by-{matrix.NumColumns} matrix");
			}
		}

		public static void CheckMultiplicationDimensions(IIndexable2D leftMatrix, IIndexable2D rightMatrix)
		{
			if (leftMatrix.NumColumns != rightMatrix.NumRows)
			{
				string message = $"Left matrix has dimensions ({leftMatrix.NumRows}x{leftMatrix.NumColumns})," +
					$" while right matrix has dimensions ({rightMatrix.NumRows}x{rightMatrix.NumColumns})";
				throw new NonMatchingDimensionsException(message);
			}
		}

		//This might be better for the methods that support transposition of the left or right matrix
		public static void CheckMultiplicationDimensions(int leftColumns, int rightRows)
		{
			if (leftColumns != rightRows)
			{
				string message = $"Left matrix/vector has {leftColumns} columns, while right matrix/vector has {rightRows} rows";
				throw new NonMatchingDimensionsException(message);
			}
		}

		public static void CheckMultiplicationDimensionsSection(IIndexable2D matrixLeft, IReadOnlyVector vectorRight,
			int vectorStart, IReadOnlyVector result, int resultStart)
		{
			if (vectorStart + matrixLeft.NumColumns > vectorRight.Length)
			{
				throw new NonMatchingDimensionsException($"The multiplied vector's length = {vectorRight.Length} must be at" +
					$" least as large as the start index = {vectorStart} + the matrix' columns = {matrixLeft.NumColumns}");
			}

			if (vectorStart + matrixLeft.NumRows > result.Length)
			{
				throw new NonMatchingDimensionsException($"The result vector's length = {result.Length} must be at least as" +
					$" large as the start index = {resultStart} + the matrix' rows = {matrixLeft.NumRows}");
			}
		}

		public static void CheckMultiplicationDimensions(IIndexable2D matrix, IReadOnlyVector lhsVector, int lhsOffset,
			IReadOnlyVector rhsVector, int rhsOffset, bool transposeMatrix)
		{
			int m, n;
			if (transposeMatrix)
			{
				m = matrix.NumColumns;
				n = matrix.NumRows;
			}
			else
			{
				m = matrix.NumRows;
				n = matrix.NumColumns;
			}

			if (lhsVector.Length < lhsOffset + n)
			{
				throw new NonMatchingDimensionsException("The left hand side vector's length must be at least as large as the" +
				  " offset plus the number of columns of the matrix (or its transpose).");
			}

			if (rhsVector.Length < rhsOffset + m)
			{
				throw new NonMatchingDimensionsException("The left hand side vector's length must be at least as large as the" +
					" offset plus the number of rows of the matrix (or its transpose).");
			}
		}

		public static void CheckMultiplicationDimensionsMatrixVector(IIndexable2D matrix, IIndexable1D lhsVector)
		{
			if (matrix.NumColumns != lhsVector.Length)
			{
				throw new NonMatchingDimensionsException(
					$"Cannot multiply ({matrix.NumRows} x {matrix.NumColumns}) * ({lhsVector.Length} x 1)");
			}
		}

		public static void CheckMultiplicationDimensionsMatrixVector(
			IIndexable2D matrix, IIndexable1D lhsVector, IIndexable1D rhsVector)
		{
			CheckMultiplicationDimensionsMatrixVector(matrix, lhsVector);
			if (matrix.NumRows != rhsVector.Length)
			{
				throw new NonMatchingDimensionsException(
					$"The result of ({matrix.NumRows} x {matrix.NumColumns}) * ({lhsVector.Length} x 1)" +
					$" must be ({matrix.NumRows} x 1), but a ({rhsVector.Length} x 1) vector was provided");
			}
		}

		public static void CheckMultiplicationDimensionsMatrixVectorTranspose(IIndexable2D matrix, IIndexable1D lhsVector)
		{
			if (matrix.NumRows != lhsVector.Length)
			{
				throw new NonMatchingDimensionsException(
					$"Cannot multiply ({matrix.NumColumns} x {matrix.NumRows}) * ({lhsVector.Length} x 1)");
			}
		}

		public static void CheckMultiplicationDimensionsMatrixVectorTranspose(
			IIndexable2D matrix, IIndexable1D lhsVector, IIndexable1D rhsVector)
		{
			CheckMultiplicationDimensionsMatrixVectorTranspose(matrix, lhsVector);
			if (matrix.NumColumns != rhsVector.Length)
			{
				throw new NonMatchingDimensionsException(
					$"The result of ({matrix.NumColumns} x {matrix.NumRows}) * ({lhsVector.Length} x 1)" +
					$" must be ({matrix.NumColumns} x 1), but a ({rhsVector.Length} x 1) vector was provided.");
			}
		}

		public static void CheckSameColDimension(IIndexable2D matrix, IReadOnlyVector vector)
		{
			if (matrix.NumColumns != vector.Length)
			{
				string message = $"Matrix has dimensions ({matrix.NumRows}x{matrix.NumColumns})," +
					$" while vector has dimensions ({vector.Length}x{matrix.NumColumns})";
				throw new NonMatchingDimensionsException(message);
			}
		}

		public static void CheckSameColDimension(IIndexable2D matrix1, IIndexable2D matrix2)
		{
			if (matrix1.NumColumns != matrix2.NumColumns)
			{
				string message = $"Matrix1 has dimensions ({matrix1.NumRows}x{matrix1.NumColumns})," +
					$" while matrix2 has dimensions ({matrix2.NumRows}x{matrix2.NumColumns})";
				throw new NonMatchingDimensionsException(message);
			}
		}

		public static void CheckSameMatrixDimensions(IIndexable2D matrix1, IIndexable2D matrix2)
		{
			if ((matrix1.NumRows != matrix2.NumRows) || (matrix1.NumColumns != matrix2.NumColumns))
			{
				string message = $"Matrix1 has dimensions ({matrix1.NumRows}x{matrix1.NumColumns})," +
					$" while matrix2 has dimensions ({matrix2.NumRows}x{matrix2.NumColumns})";
				throw new NonMatchingDimensionsException(message);
			}
		}

		public static void CheckSameRowDimension(IIndexable2D matrix, IReadOnlyVector vector)
		{
			if (matrix.NumRows != vector.Length)
			{
				string message = $"Matrix has dimensions ({matrix.NumRows}x{matrix.NumColumns})," +
					$" while vector has dimensions ({vector.Length}x{matrix.NumColumns})";
				throw new NonMatchingDimensionsException(message);
			}
		}

		public static void CheckSameRowDimension(IIndexable2D matrix1, IIndexable2D matrix2)
		{
			if (matrix1.NumRows != matrix2.NumRows)
			{
				string message = $"Matrix1 has dimensions ({matrix1.NumRows}x{matrix1.NumColumns})," +
					$" while matrix2 has dimensions ({matrix2.NumRows}x{matrix2.NumColumns})";
				throw new NonMatchingDimensionsException(message);
			}
		}

		public static void CheckSquare(IIndexable2D matrix)
		{
			if (matrix.NumRows != matrix.NumColumns)
			{
				throw new NonMatchingDimensionsException(
					$"The matrix must be square, but was {matrix.NumRows}-by-{matrix.NumColumns}");
			}
		}

		public static void CheckSquare(int numRows, int numColumns)
		{
			if (numRows != numColumns)
			{
				throw new NonMatchingDimensionsException($"The matrix must be square, but was {numRows}-by-{numColumns}");
			}
		}

		public static void CheckSquareLinearSystemDimensions(IIndexable2D matrix, IIndexable1D lhsVector, IIndexable1D rhsVector)
		{
			CheckSquareLinearSystemDimensions(matrix.NumRows, matrix.NumColumns, lhsVector.Length, rhsVector.Length);
		}

		public static void CheckSquareLinearSystemDimensions(int numMatrixRows, int numMatrixColumns, int lhsLength, int rhsLength)
		{
			bool compatible = numMatrixRows == numMatrixColumns;
			compatible &= numMatrixColumns == lhsLength;
			compatible &= numMatrixRows == rhsLength;
			if (!compatible)
			{
				throw new NonMatchingDimensionsException(
					$"The matrix rows ({numMatrixRows}), matrix columns ({numMatrixColumns})," +
					$" left-hand-side vector length ({lhsLength}) and" +
					$" right-hand-side vector length ({rhsLength}) must be the same.");
			}
		}

		public static void CheckSubvectorDimensions(IIndexable1D vector, int startIndex, int subvectorLength)
		{
			if (startIndex + subvectorLength > vector.Length)
			{
				throw new NonMatchingDimensionsException("The entries to access exceed the vector's length");
			}
		}

		public static void CheckSystemSolutionDimensions(IIndexable2D matrix, IReadOnlyVector rhsVector)
		{
			if (matrix.NumRows != rhsVector.Length)
			{
				string message = $"Matrix has dimensions ({matrix.NumRows}x{matrix.NumColumns})," +
					$" while the right hand side vector has dimensions ({rhsVector.Length}x1)";
				throw new NonMatchingDimensionsException(message);
			}
		}

		public static void CheckSystemSolutionDimensions(int matrixRows, int rhsVectorLength) //TODO: this should be deleted
		{
			if (matrixRows != rhsVectorLength)
			{
				string message = $"Matrix has {matrixRows} rows, while the right hand side vector has"
					+ $" dimensions ({rhsVectorLength}x1)";
				throw new NonMatchingDimensionsException(message);
			}
		}

		public static void CheckVectorDimensions(int vector1Length, int vector2Length)
		{
			if (vector1Length != vector2Length)
			{
				string message = $"Vector1 has length of {vector1Length}, while vector2 has length of {vector2Length}";
				throw new NonMatchingDimensionsException(message);
			}
		}

		public static void CheckVectorDimensions(double[] vector1, double[] vector2)
		{
			if (vector1.Length != vector2.Length)
			{
				string message = $"Vector1 has length of {vector1.Length}, while vector2 has length of {vector2.Length}";
				throw new NonMatchingDimensionsException(message);
			}

		}

		public static void CheckVectorDimensions(IIndexable1D vector1, IIndexable1D vector2)
		{
			if (vector1.Length != vector2.Length)
			{
				string message = $"Vector1 has length of {vector1.Length}, while vector2 has length of {vector2.Length}";
				throw new NonMatchingDimensionsException(message);
			}
		}
	}
}
