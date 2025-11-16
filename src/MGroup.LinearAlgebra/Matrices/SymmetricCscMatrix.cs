//TODO: Should this be removed? I need it to provide analyzers with a matrix.
namespace MGroup.LinearAlgebra.Matrices
{
	using System;

	using MGroup.LinearAlgebra.Commons;
	using MGroup.LinearAlgebra.Exceptions;
	using MGroup.LinearAlgebra.Implementations;
	using MGroup.LinearAlgebra.Matrices.Builders;
	using MGroup.LinearAlgebra.Reduction;
	using MGroup.LinearAlgebra.Vectors;

	/// <summary>
	/// Symmetric sparse matrix in Compressed Sparse Columns format, with only the non-zero entries of the upper triangle being
	/// explicitly stored in column major order. This matrix format is best used for factorizations using SuiteSparse
	/// and CSparse libraries.
	/// </summary>
	public sealed class SymmetricCscMatrix : ValuesBackedMatrix<SymmetricCscMatrix>
	{
		/// <summary>
		/// values.Length = number of non zeros in upper triangle.
		/// </summary>
		private readonly double[] values;

		/// <summary>
		/// rowIndices.Length = number of non zeros in upper triangle.
		/// </summary>
		private readonly int[] rowIndices;

		/// <summary>
		/// colOffsets.Length = number of rows/columns +1.
		/// colOffsets[colOffsets.Length-1] = number of non zeros in upper triangle.
		/// </summary>
		private readonly int[] colOffsets;

		private SymmetricCscMatrix(int order, double[] values, int[] rowIndices, int[] colOffsets)
		{
			this.values = values;
			this.NumColumns = order;
			this.NumRows = order;
			this.NumNonZerosUpper = values.Length;
			this.rowIndices = rowIndices;
			this.colOffsets = colOffsets;
		}

		public override MatrixSymmetry MatrixSymmetry => MatrixSymmetry.Symmetric;

		public override int NumColumns { get; }

		/// <summary>
		/// The number of the upper triangle's non zero entries. These are the only ones being explicitly stored.
		/// </summary>
		public int NumNonZerosUpper { get; }

		public override int NumRows { get; }

		/// <summary>
		/// The internal array that stores the index into the arrays <see cref="ValuesBackedMatrix{TMatrix}.RawValues"/> and
		/// <see cref="RawRowIndices"/> of the first entry of each column. Its length is equal to
		/// <paramref name="NumColumns"/> + 1.
		/// The last entry is the number of the upper triangle's non-zero entries, which must be equal to
		/// <see cref="ValuesBackedMatrix{TMatrix}.RawValues"/>.Length == <see cref="RawRowIndices"/>.Length.
		/// It should only be used for passing the raw array to linear algebra libraries.
		/// </summary>
		public int[] RawColOffsets => colOffsets;

		/// <summary>
		/// The internal array that stores the row indices of the non-zero entries in <see cref="RawValues"/>.
		/// Its length is equal to the number of the upper triangle's non-zero entries.
		/// It should only be used for passing the raw array to linear algebra libraries.
		/// </summary>
		public int[] RawRowIndices => rowIndices;

		public override double[] RawValues => values;

		public override double this[int rowIdx, int colIdx]
		{
			get
			{
				int offset = FindOffsetOf(rowIdx, colIdx);
				if (offset >= 0) return values[offset];
				else return 0.0;
			}

			set
			{
				int offset = FindOffsetOf(rowIdx, colIdx);
				if (offset >= 0) values[offset] = value;
				else throw new SparsityPatternModifiedException($"Cannot write to zero entry ({rowIdx}, {colIdx}).");
			}
		}

		/// <summary>
		/// Initializes a new <see cref="SymmetricCscMatrix"/> with the specified dimensions and the provided arrays
		/// (<paramref name="values"/>, <paramref name="rowIndices"/> and <paramref name="colOffsets"/>) as its internal data.
		/// </summary>
		/// <param name="order">The number of rows /columns of the new matrix.</param>
		/// <param name="values">
		/// Array that contains the non-zero entries of the upper triangle. It must have the same length as
		/// <paramref name="rowIndices"/>. The non-zero entries of each column must appear consecutively in
		/// <paramref name="values"/>. They can also be sorted in increasing order of their row indices, which speeds up
		/// subsequent operations.
		/// </param>
		/// <param name="rowIndices">
		/// Array that contains the row indices of the upper triangle's non-zero entries. It must have the same length as
		/// <paramref name="values"/>. There is an 1 to 1 matching between these two arrays: <paramref name="rowIndices"/>[i]
		/// is the row index of the entry <paramref name="values"/>[i]. Also:
		/// 0 &lt;= <paramref name="rowIndices"/>[i] &lt; <paramref name="numRows"/>.
		/// </param>
		/// <param name="colOffsets">
		/// Array that contains the index of the first entry of each column into the arrays <paramref name="values"/> and
		/// <paramref name="rowIndices"/>. Its length is <paramref name="numRows"/> + 1. The last entry is the number of
		/// non-zero entries, which must be equal to the length of <paramref name="values"/>and <paramref name="rowIndices"/>.
		/// </param>
		/// <param name="checkInput">
		/// If true, the provided arrays will be checked to make sure they are valid symmetric CSC arrays, which is safer.
		/// If false, no such check will take place, which is faster.
		/// </param>
		public static SymmetricCscMatrix CreateFromArrays(int order, double[] values, int[] rowIndices, int[] colOffsets,
			bool checkInput)
		{
			int nnz = colOffsets[colOffsets.Length - 1];
			if (checkInput)
			{
				if (colOffsets.Length != order + 1)
				{
					throw new ArgumentException("The length of the symmetric CSC column offsets array must be equal to the order"
						+ " of the matrix + 1, but was " + colOffsets.Length);
				}
				if ((nnz != values.Length) || (nnz != rowIndices.Length))
				{
					throw new ArgumentException("Mismatch in dimensions of the symmetric CSC arrays. Check that"
						+ " colOffsets.Length = number of rows/columns + 1 and that colOffsets[colOffsets.Length-1]"
						+ " = values.Length = rowIndices.Length = number of non zeros in upper triangle");
				}
			}
			return new SymmetricCscMatrix(order, values, rowIndices, colOffsets);
		}

		public override SymmetricCscMatrix CopyAsSameType(bool copyIndexingArrays)
		{
			var valuesCopy = new double[values.Length];
			Array.Copy(values, valuesCopy, values.Length);
			if (!copyIndexingArrays)
			{
				return new SymmetricCscMatrix(NumColumns, valuesCopy, rowIndices, colOffsets);
			}
			else
			{
				var rowIndicesCopy = new int[rowIndices.Length];
				Array.Copy(rowIndices, rowIndicesCopy, rowIndices.Length);
				var colOffsetsCopy = new int[colOffsets.Length];
				Array.Copy(colOffsets, colOffsetsCopy, colOffsets.Length);
				return new SymmetricCscMatrix(NumColumns, valuesCopy, rowIndicesCopy, colOffsetsCopy);
			}
		}

		/// <summary>
		/// Creates a CSR matrix with the same entries as this. The returned CSR matrix is symmetric, but stores both triangles.
		/// </summary>
		/// <returns>A CSR matrix with the same entries as this.</returns>
		public CsrMatrix ConvertToCsr()
		{
			var dok = DokRowMajor.CreateEmpty(NumColumns, NumColumns);
			for (int j = 0; j < NumColumns; ++j)
			{
				int colStart = colOffsets[j];
				int colEnd = colOffsets[j + 1];
				for (int k = colStart; k < colEnd; ++k)
				{
					int i = rowIndices[k];
					double val = values[k];
					dok[i, j] = val;

					if (i != j)
					{
						dok[j, i] = val;
					}
				}
			}

			return dok.BuildCsrMatrix(true);
		}

		public override Matrix CopyToFullMatrix()
		{
			Matrix fullMatrix = Matrix.CreateZero(this.NumRows, this.NumColumns);
			for (int j = 0; j < this.NumColumns; ++j) //Column major order
			{
				int colCurrent = colOffsets[j];
				int colNext = colOffsets[j + 1]; // TODO: 1 of the two accesses can be removed
				for (int k = colCurrent; k < colNext; ++k)
				{
					int i = rowIndices[k];
					double val = values[k];
					fullMatrix[i, j] = val;
					fullMatrix[j, i] = val;
				}
			}
			return fullMatrix;
		}

		public override SymmetricCscMatrix CreateZeroMatrixSame()
		{
			var resultValues = new double[values.Length];
			return new SymmetricCscMatrix(NumColumns, resultValues, rowIndices, colOffsets);
		}

		public override Vector GetRow(int rowIndex) => GetColumn(rowIndex);

		public override bool HasSameFormat(SymmetricCscMatrix other)
		{
			return (this.rowIndices == other.rowIndices) && (this.colOffsets == other.colOffsets);
		}

		public override IVector Multiply(IReadOnlyVector vector, bool transposeThis = false)
		{
			var result = Vector.CreateZero(NumRows);
			CsrMultiplications.SymmetricCsrTimesVector(NumRows, values, colOffsets, rowIndices, vector, result.RawData);
			return result;
		}

		public override void MultiplyIntoResult(IReadOnlyVector lhsVector, IVector rhsVector, bool transposeThis)
		{
			rhsVector.Clear(); // TODO: add this as an optional flag or convert the method to axpy like.
			if (rhsVector is Vector denseVector)
			{
				CsrMultiplications.SymmetricCsrTimesVector(
					NumRows, values, colOffsets, rowIndices, lhsVector, denseVector.RawData);
			}
			else
			{
				var temp = Vector.CreateZero(NumRows);
				CsrMultiplications.SymmetricCsrTimesVector(NumRows, values, colOffsets, rowIndices, lhsVector, temp.RawData);
				rhsVector.CopyFrom(temp);
			}
		}

		/// <summary>
		/// Matrix-vector multiplication, with the vector on the right: matrix * vector or transpose(matrix) * vector.
		/// </summary>
		/// <param name="vector">A vector with length equal to <see cref="NumColumns"/>.</param>
		/// <param name="transposeThis">Set to true to transpose this (the left matrix). Unless the transpose matrix is used in
		///     more than one multiplications, setting this flag to true is usually preferable to creating the transpose.</param>
		/// <returns></returns>
		public Vector MultiplyRight(Vector vector, bool transposeThis = false)
		{
			var result = Vector.CreateZero(NumRows);
			CsrMultiplications.SymmetricCsrTimesVector(NumRows, values, colOffsets, rowIndices, vector.RawData, result.RawData);
			return result;
		}

		public override double Reduce(double identityValue, ProcessEntry processEntry, ProcessZeros processZeros, Finalize finalize)
		{
			double aggregator = identityValue;
			int numNonZeros = 0;
			for (int j = 0; j < NumColumns; ++j)
			{
				int colStart = colOffsets[j]; //inclusive
				int colEnd = colOffsets[j + 1]; //exclusive
				for (int k = colStart; k < colEnd; ++k)
				{

					if (rowIndices[k] == j)
					{
						aggregator = processEntry(values[k], aggregator);
						++numNonZeros;
					}
					else // Do the above twice for entries not on the diagonal
					{
						aggregator = processEntry(values[k], aggregator);
						aggregator = processEntry(values[k], aggregator);
						numNonZeros += 2;
					}
				}
			}
			aggregator = processZeros(NumRows * NumColumns - numNonZeros, aggregator);
			return finalize(aggregator);
		}

		public override IMatrix Transpose() => CopyAsSameType(true);

		/// <summary>
		/// For an entry (i,j), returns the offset into values and rowIndices arrays or -1 of the entry corresponds to a
		/// structural zero.
		/// </summary>
		/// <param name="rowIdx"></param>
		/// <param name="colIdx"></param>
		private int FindOffsetOf(int rowIdx, int colIdx)
		{
			if (rowIdx > colIdx)
			{
				int swap = rowIdx;
				rowIdx = colIdx;
				colIdx = swap;
			}
			int colStart = colOffsets[colIdx];
			int colEnd = colOffsets[colIdx + 1];
			for (int k = colStart; k < colEnd; ++k) //Only scan the nnz entries for the given column
			{
				if (rowIndices[k] == rowIdx) return k;
			}
			return -1;
		}
	}
}
