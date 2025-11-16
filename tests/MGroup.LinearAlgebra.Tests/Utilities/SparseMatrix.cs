namespace MGroup.LinearAlgebra.Tests.Utilities
{
	using System;
	using System.Collections.Generic;
	using System.Diagnostics;
	using System.Linq;
	using System.Text;
	using System.Threading.Tasks;

	using MGroup.LinearAlgebra.Exceptions;
	using MGroup.LinearAlgebra.Matrices;
	using MGroup.LinearAlgebra.Matrices.Builders;
	using MGroup.LinearAlgebra.Vectors;

	public class SparseMatrix
	{
		private Dictionary<int, double>[] rows;

		public SparseMatrix(int numRows, int numColumns)
		{
			NumRows = numRows;
			NumColumns = numColumns;
			rows = new Dictionary<int, double>[numRows];
			for (int i = 0; i < numRows; ++i) rows[i] = new Dictionary<int, double>();
		}

		public int NumRows { get; }

		public int NumColumns { get; }

		public double this[int rowIdx, int colIdx]
		{
			get
			{
				CheckIndex(rowIdx, colIdx);
				if (rows[rowIdx].TryGetValue(colIdx, out double val)) return val;
				else return 0.0;
			}
			set
			{
				CheckIndex(rowIdx, colIdx);
				rows[rowIdx][colIdx] = value;
			}
		}

		public static SparseMatrix CreateFromMatrix(IIndexable2D original, double dropTolerance = 1E-7)
		{
			var result = new SparseMatrix(original.NumRows, original.NumColumns);
			for (int i = 0; i < original.NumRows; ++i)
			{
				for (int j = 0; j < original.NumColumns; ++j)
				{
					double val = original[i, j];
					if (val > dropTolerance)
					{
						result[i, j] = val;
					}
				}
			}
			return result;
		}

		public static SparseMatrix CreateFromMatrix(ISparseMatrix original)
		{
			var result = new SparseMatrix(original.NumRows, original.NumColumns);
			foreach ((int row, int col, double val) in original.EnumerateNonZeros())
			{
				result[row, col] = val;
			}
			return result;
		}

		public static SparseMatrix operator +(SparseMatrix matrix1, SparseMatrix matrix2)
			=> matrix1.LinearCombination(1.0, matrix2, 1.0);

		public static SparseMatrix operator -(SparseMatrix matrix1, SparseMatrix matrix2)
			=> matrix1.LinearCombination(1.0, matrix2, -1.0);

		public static SparseMatrix operator *(double scalar, SparseMatrix matrix) => matrix.Scale(scalar);

		public static SparseMatrix operator *(SparseMatrix matrix, double scalar) => matrix.Scale(scalar);

		public static Vector operator *(SparseMatrix matrix, Vector vector) => matrix.Multiply(vector);

		public void AddSubmatrix(IIndexable2D subMatrix, int[] subMatrixRows, int[] globalMatrixRows,
			int[] subMatrixCols, int[] globalMatrixCols)
		{
			int numRows = subMatrixRows.Length;
			int numCols = subMatrixCols.Length;
			Debug.Assert(numRows == globalMatrixRows.Length);
			Debug.Assert(numCols == globalMatrixCols.Length);

			for (int i = 0; i < numRows; ++i)
			{
				int subRow = subMatrixRows[i];
				int globalRow = globalMatrixRows[i];
				Debug.Assert((globalRow >= 0) && (globalRow < NumRows));

				for (int j = 0; j < numCols; ++j)
				{
					int subCol = subMatrixCols[j];
					int globalCol = globalMatrixCols[j];
					Debug.Assert((globalCol >= 0) && (globalCol < NumColumns));

					double subVal = subMatrix[subRow, subCol];
					rows[globalRow].TryGetValue(globalCol, out double oldGlobalVal); // default value = 0.0, if the entry is new
					rows[globalRow][globalCol] = subVal + oldGlobalVal;
				}
			}
		}

		public SparseMatrix Copy()
		{
			var clone = new SparseMatrix(NumRows, NumColumns);
			for (int i = 0; i < NumRows; ++i)
			{
				Dictionary<int, double> wholeRow = rows[i];
				clone.rows[i] = new Dictionary<int, double>(this.rows[i]);
			}

			return clone;
		}

		public IEnumerable<(int row, int col, double value)> EnumerateNonZeros()
		{
			for (int i = 0; i < NumRows; ++i)
			{
				foreach (var colVal in rows[i])
				{
					yield return (i, colVal.Key, colVal.Value);
				}
			}
		}

		public SparseMatrix ExtractDiagonal()
		{
			var result = new SparseMatrix(NumRows, NumColumns);
			for (int i = 0; i < NumRows; ++i)
			{
				Dictionary<int, double> wholeRow = rows[i];
				double val = 0.0;
				wholeRow.TryGetValue(i, out val);
				result[i, i] = val;
			}
			return result;
		}

		public SparseMatrix ExtractLowerTriangle()
		{
			var result = new SparseMatrix(NumRows, NumColumns);
			for (int i = 0; i < NumRows; ++i)
			{
				Dictionary<int, double> wholeRow = rows[i];
				foreach (KeyValuePair<int, double> colValuePair in wholeRow)
				{
					int j = colValuePair.Key;
					if (j < i)
					{
						double val = colValuePair.Value;
						result[i, j] = val;
					}
				}
			}
			return result;
		}

		public SparseMatrix ExtractLowerTriangleAndDiagonal()
		{
			var result = new SparseMatrix(NumRows, NumColumns);
			for (int i = 0; i < NumRows; ++i)
			{
				Dictionary<int, double> wholeRow = rows[i];
				foreach (KeyValuePair<int, double> colValuePair in wholeRow)
				{
					int j = colValuePair.Key;
					if (j <= i)
					{
						double val = colValuePair.Value;
						result[i, j] = val;
					}
				}
			}
			return result;
		}

		public SparseMatrix ExtractSubmatrix(int[] rowsToKeep, int[] colsToKeep)
		{
			var oldToNewRows = new Dictionary<int, int>();
			for (int i = 0; i < rowsToKeep.Length; ++i)
			{
				oldToNewRows[rowsToKeep[i]] = i;
			}

			var oldToNewCols = new Dictionary<int, int>();
			for (int j = 0; j < colsToKeep.Length; ++j)
			{
				oldToNewCols[colsToKeep[j]] = j;
			}

			var result = new SparseMatrix(rowsToKeep.Length, colsToKeep.Length);
			for (int I = 0; I < this.NumRows; ++I) // Traverse the existing sparse matrix and copy only the requested entries
			{
				bool keepRow = oldToNewRows.TryGetValue(I, out int i);
				if (!keepRow)
				{
					continue;
				}

				foreach (var colValPair in this.rows[I])
				{
					int J = colValPair.Key;
					bool keepCol = oldToNewCols.TryGetValue(J, out int j);
					if (!keepCol)
					{
						continue;
					}

					double val = colValPair.Value;
					result[i, j] = val;
				}
			}

			return result;
		}

		public SparseMatrix ExtractUpperTriangle()
		{
			var result = new SparseMatrix(NumRows, NumColumns);
			for (int i = 0; i < NumRows; ++i)
			{
				Dictionary<int, double> wholeRow = rows[i];
				foreach (KeyValuePair<int, double> colValuePair in wholeRow)
				{
					int j = colValuePair.Key;
					if (j > i)
					{
						double val = colValuePair.Value;
						result[i, j] = val;
					}
				}
			}
			return result;
		}

		public SparseMatrix ExtractUpperTriangleAndDiagonal()
		{
			var result = new SparseMatrix(NumRows, NumColumns);
			for (int i = 0; i < NumRows; ++i)
			{
				Dictionary<int, double> wholeRow = rows[i];
				foreach (KeyValuePair<int, double> colValuePair in wholeRow)
				{
					int j = colValuePair.Key;
					if (j >= i)
					{
						double val = colValuePair.Value;
						result[i, j] = val;
					}
				}
			}
			return result;
		}

		public SparseMatrix LinearCombination(double thisCoeff, SparseMatrix other, double otherCoeff)
		{
			CheckSameDimensions(other);
			SparseMatrix result = this.Copy();
			for (int i = 0; i < NumRows; ++i)
			{
				Dictionary<int, double> otherRow = other.rows[i];
				Dictionary<int, double> resultRow = result.rows[i];
				foreach (KeyValuePair<int, double> colValuePair in otherRow)
				{
					int j = colValuePair.Key;
					double otherVal = colValuePair.Value;
					resultRow.TryGetValue(j, out double resultVal);
					resultRow[j] = thisCoeff * resultVal + otherCoeff * otherVal;
				}
			}

			return result;
		}

		public Vector Multiply(Vector vector)
		{
			CheckVectorMult(vector);
			var result = new double[this.NumRows];
			for (int i = 0; i < NumRows; ++i)
			{
				double dot = 0.0;
				foreach (KeyValuePair<int, double> colValPair in rows[i])
				{
					int j = colValPair.Key;
					double val = colValPair.Value;
					dot += val * vector[j];
				}
				result[i] = dot;
			}
			return Vector.CreateFromArray(result, false);
		}

		public Vector MultiplyTranspose(Vector vector)
		{
			CheckSystemSolution(vector);
			var result = new double[this.NumColumns];
			for (int i = 0; i < NumRows; ++i)
			{
				double scale = vector[i];
				foreach (KeyValuePair<int, double> colValPair in rows[i])
				{
					int j = colValPair.Key;
					double val = colValPair.Value;
					result[j] += scale * val;
				}
			}
			return Vector.CreateFromArray(result, false);
		}

		public SparseMatrix Scale(double scalar)
		{
			var result = new SparseMatrix(NumRows, NumColumns);
			for (int i = 0; i < NumRows; ++i)
			{
				foreach (KeyValuePair<int, double> colValPair in rows[i])
				{
					int j = colValPair.Key;
					double val = colValPair.Value;
					result[i, j] = scalar * val;
				}
			}
			return result;
		}

		public Vector SolveDiagonal(Vector rhs)
		{
			CheckSystemSolution(rhs);
			int n = rhs.Length;
			var solution = new double[n];
			for (int i = 0; i < rhs.Length; ++i)
			{
				solution[i] = rhs[i] / this[i, i];
			}
			return Vector.CreateFromArray(solution);
		}

		public Vector SolveForwardSubstitution(Vector rhs)
		{
			CheckSystemSolution(rhs);
			int n = rhs.Length;
			var solution = new double[n];
			for (int i = 0; i < rhs.Length; ++i)
			{
				double denom = 0.0;
				double sum = 0.0;
				foreach (KeyValuePair<int, double> colValPair in rows[i])
				{
					int j = colValPair.Key;
					double val = colValPair.Value;
					if (j == i)
					{
						denom = val;
					}
					else if (j < i)
					{
						sum += val * solution[j];
					}
					else
					{
						throw new Exception("This matrix is not lower triangular");
					}
				}
				solution[i] = (rhs[i] - sum) / denom;
			}
			return Vector.CreateFromArray(solution);
		}

		public Vector SolveBackSubstitution(Vector rhs)
		{
			CheckSystemSolution(rhs);
			int n = rhs.Length;
			var solution = new double[n];
			for (int i = rhs.Length - 1; i >= 0; --i)
			{
				double denom = 0.0;
				double sum = 0.0;
				foreach (KeyValuePair<int, double> colValPair in rows[i])
				{
					int j = colValPair.Key;
					double val = colValPair.Value;
					if (j == i)
					{
						denom = val;
					}
					else if (j > i)
					{
						sum += val * solution[j];
					}
					else
					{
						throw new Exception("This matrix is not upper triangular");
					}
				}
				solution[i] = (rhs[i] - sum) / denom;
			}
			return Vector.CreateFromArray(solution);
		}

		private void CheckIndex(int row, int col)
		{
			if (row < 0 || row >= NumRows)
			{
				throw new IndexOutOfRangeException($"Row index must be in [0, {NumRows}], but was {row}.");
			}

			if (col < 0 || col >= NumColumns)
			{
				throw new IndexOutOfRangeException($"Colum index must be in [0, {NumColumns}], but was {col}.");
			}
		}

		private void CheckSameDimensions(SparseMatrix other)
		{
			if ((other.NumRows != this.NumRows) || (other.NumColumns != this.NumColumns))
			{
				throw new NonMatchingDimensionsException($"The other matrix ({other.NumRows}x{other.NumColumns}) must have the " +
					$"same dimensions as this matrix ({this.NumRows}x{this.NumColumns}).");
			}
		}

		private void CheckSystemSolution(IReadOnlyVector vector)
		{
			if (NumRows != NumColumns)
			{
				throw new NonMatchingDimensionsException(
					$"The matrix must be square but has {NumRows} rows and {NumColumns} columns.");
			}
			if (vector.Length != NumRows)
			{
				throw new NonMatchingDimensionsException($"The vector must have {NumColumns} entries, but has {vector.Length}.");
			}
		}

		private void CheckVectorMult(IReadOnlyVector vector)
		{
			if (vector.Length != NumColumns)
			{
				throw new NonMatchingDimensionsException($"The vector must have {NumColumns} entries, but has {vector.Length}.");
			}
		}
	}
}
