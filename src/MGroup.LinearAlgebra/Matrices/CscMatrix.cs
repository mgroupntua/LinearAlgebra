//TODO: try to make general versions of row major and col major multiplication. The lhs matrix/vector will be supplied by the 
// caller, depending on if it is a matrix, a transposed matrix, a vector, etc. Compare its performance with the verbose code and 
// also try inlining. C preprocessor macros would be actually useful here. Otherwise, move all that boilerplate code to a 
// CSRStrategies static class.
//TODO: In matrix-matrix/vector multiplications: perhaps I should work with a column major array directly instead of an output   
//      Matrix and an array instead of an output Vector.
//TODO: perhaps optimizations if (other is Matrix) are needed, to directly index into its raw col major array.
//      The access paterns are always the same.
//TODO: The implementations of this class should call transposed operations on a backing CSR matrix.
namespace MGroup.LinearAlgebra.Matrices
{
	using System;
	using System.Collections.Generic;

	using MGroup.LinearAlgebra.Commons;
	using MGroup.LinearAlgebra.Exceptions;
	using MGroup.LinearAlgebra.Output.Formatting;
	using MGroup.LinearAlgebra.Reduction;
	using MGroup.LinearAlgebra.Vectors;

	using static MGroup.LinearAlgebra.LibrarySettings;

	/// <summary>
	/// Sparse matrix stored in Compressed Sparse Columns format (3-array version). The CSR format is optimized for matrix-vector 
	/// and matrix-matrix multiplications, where the CSC matrix is on the left transposed or on the right untransposed. The other
	/// multiplicationss are more efficient using <see cref="CsrMatrix"/>. To build a <see cref="CscMatrix"/> conveniently, 
	/// use <see cref="Builders.DokColMajor"/>.
	/// </summary>
	public sealed class CscMatrix : ValuesBackedMatrix<CscMatrix>, ISparseMatrix
	{
		private const int zeroEntryOffset = -1;

		private readonly double[] values;
		private readonly int[] rowIndices;
		private readonly int[] colOffsets;

		private CscMatrix(int numRows, int numCols, double[] values, int[] rowIndices, int[] colOffsets)
		{
			this.values = values;
			this.rowIndices = rowIndices;
			this.colOffsets = colOffsets;
			this.NumRows = numRows;
			this.NumColumns = numCols;
		}

		public override int NumColumns { get; }

		/// <summary>
		/// The number of non zero entries of the matrix.
		/// </summary>
		public int NumNonZeros => rowIndices.Length;

		public override int NumRows { get; }

		/// <summary>
		/// The internal array that stores the index into the arrays <see cref="RawValues"/> and <see cref="RawRowIndices"/> of
		/// the first entry of each column. Its length is equal to <paramref name="NumColumns"/> + 1.
		/// The last entry is the number of non-zero entries, which must be equal to
		/// <see cref="ValuesBackedMatrix{TMatrix}.RawValues"/>.Length == <see cref="RawRowIndices"/>.Length.
		/// It should only be used for passing the raw array to linear algebra libraries.
		/// </summary>
		public int[] RawColOffsets => colOffsets;

		/// <summary>
		/// The internal array that stores the row indices of the non-zero entries in
		/// <see cref="ValuesBackedMatrix{TMatrix}.RawValues"/>.
		/// Its length is equal to the number of non-zero entries.
		/// It should only be used for passing the raw array to linear algebra libraries.
		/// </summary>
		public int[] RawRowIndices => rowIndices;

		public override double[] RawValues => values;


		public override double this[int rowIdx, int colIdx]
		{
			get
			{
				int entryOffset = FindOffsetOf(rowIdx, colIdx);
				if (entryOffset == zeroEntryOffset) return 0.0;
				else return values[entryOffset];
			}

			set
			{
				int entryOfsset = FindOffsetOf(rowIdx, colIdx);
				if (entryOfsset == zeroEntryOffset) throw new SparsityPatternModifiedException(
					$"Cannot write to zero entry ({rowIdx}, {colIdx}).");
				else values[entryOfsset] = value;
			}
		}

		/// <summary>
		/// Initializes a new <see cref="CscMatrix"/> with the specified dimensions and the provided arrays 
		/// (<paramref name="values"/>, <paramref name="rowIndices"/> and <paramref name="colOffsets"/>) as its internal data.
		/// </summary>
		/// <param name="numRows">The number of rows of the new matrix.</param>
		/// <param name="numCols">The number of columns of the new matrix.</param>
		/// <param name="values">Array that contains the non-zero entries. It must have the same length 
		///     as <paramref name="rowIndices"/>. The non-zero entries of each column must appear consecutively in 
		///     <paramref name="values"/>. They can also be sorted in increasing order of their row indices, which speeds up
		///     subsequent operations.</param>
		/// <param name="rowIndices">Array that contains the row indices of the non-zero entries. It must have the same 
		///     length as <paramref name="values"/>. There is an 1 to 1 matching between these two arrays: 
		///     <paramref name="rowIndices"/>[i] is the row index of the entry <paramref name="values"/>[i]. Also:
		///     0 &lt;= <paramref name="rowIndices"/>[i] &lt; <paramref name="numRows"/>.</param>
		/// <param name="colOffsets">Array that contains the index of the first entry of each column into the arrays 
		///     <paramref name="values"/> and <paramref name="rowIndices"/>. Its length is <paramref name="numRows"/> + 1. The 
		///     last entry is the number of non-zero entries, which must be equal to the length of <paramref name="values"/> 
		///     and <paramref name="rowIndices"/>.</param>
		/// <param name="checkInput">If true, the provided arrays will be checked to make sure they are valid CSC arrays, which 
		///     is safer. If false, no such check will take place, which is faster.</param>
		public static CscMatrix CreateFromArrays(int numRows, int numCols, double[] values, int[] rowIndices, int[] colOffsets,
			bool checkInput)
		{
			if (checkInput)
			{
				if (colOffsets.Length != numCols + 1)
				{
					throw new ArgumentException("The length of the CSC column offsets array must be equal to the number of"
						+ " rows + 1, but was " + colOffsets.Length);
				}
				if (values.Length != rowIndices.Length)
				{
					throw new ArgumentException("The length of the CSC values and row indices arrays must be equal (and equal"
						+ $" to the number of non zero entries), but were {values.Length} and {rowIndices.Length} respectively");
				}
				if (colOffsets[0] != 0)
				{
					throw new ArgumentException("The first entry of the CSC column offsets array must be 0, but was "
						+ colOffsets[0]);
				}
				if (colOffsets[colOffsets.Length - 1] != values.Length)
				{
					throw new ArgumentException("The last entry of the CSC column offsets array must be equal to the number of"
						+ " non zero entries, but was " + colOffsets[colOffsets.Length - 1]);
				}
			}
			return new CscMatrix(numRows, numCols, values, rowIndices, colOffsets);
		}

		#region operators (use extension operators when they become available)
		/// <summary>
		/// Performs the matrix-vector multiplication: result = <paramref name="vectorLeft"/> * <paramref name="matrixRight"/>.
		/// If <paramref name="matrixRight"/> is m1-by-n1 and <paramref name="vectorLeft"/> has length = n2, then m1 must be 
		/// equal to n2. The result will be a vector with length = n1, written to a new <see cref="Vector"/> instance.
		/// </summary>
		/// <param name="vectorLeft">The <see cref="Vector"/> operand on the left. It can be considered as a row vector.</param>
		/// <param name="matrixRight">The <see cref="CscMatrix"/> operand on the right.</param>
		/// <exception cref="NonMatchingDimensionsException">Thrown if <paramref name="matrixRight"/>.<see cref="NumRows"/> is 
		///     different than <paramref name="vectorLeft"/>.<see cref="Vector.Length"/>.</exception>
		public static Vector operator *(Vector vectorLeft, CscMatrix matrixRight)
			=> matrixRight.Multiply(vectorLeft, true);
		#endregion

		public override CscMatrix CopyAsSameType(bool copyIndexingData)
		{
			var valuesCopy = new double[this.values.Length];
			Array.Copy(this.values, valuesCopy, this.values.Length);

			if (!copyIndexingData) return new CscMatrix(NumRows, NumColumns, valuesCopy, this.rowIndices, this.colOffsets);
			else
			{
				var rowIndicesCopy = new int[this.rowIndices.Length];
				Array.Copy(this.rowIndices, rowIndicesCopy, this.rowIndices.Length);
				var colOffsetsCopy = new int[this.colOffsets.Length];
				Array.Copy(this.colOffsets, colOffsetsCopy, this.colOffsets.Length);
				return new CscMatrix(NumRows, NumColumns, valuesCopy, rowIndicesCopy, colOffsetsCopy);
			}
		}

		public override Matrix CopyToFullMatrix()
		{
			Matrix fullMatrix = Matrix.CreateZero(this.NumRows, this.NumColumns);
			for (int j = 0; j < this.NumColumns; ++j) //Row major order
			{
				int colStart = colOffsets[j]; //inclusive
				int colEnd = colOffsets[j + 1]; //exclusive
				for (int k = colStart; k < colEnd; ++k)
				{
					fullMatrix[rowIndices[k], j] = values[k];
				}
			}
			return fullMatrix;
		}

		public int CountNonZeros() => values.Length;

		public override CscMatrix CreateZeroMatrixSame()
		{
			var resultValues = new double[values.Length];
			return new CscMatrix(NumRows, NumColumns, resultValues, rowIndices, colOffsets);
		}

		public IEnumerable<(int row, int col, double value)> EnumerateNonZeros()
		{
			for (int j = 0; j < NumColumns; ++j)
			{
				int colStart = colOffsets[j]; //inclusive
				int colEnd = colOffsets[j + 1]; //exclusive
				for (int k = colStart; k < colEnd; ++k)
				{
					yield return (rowIndices[k], j, values[k]);
				}
			}
		}

		public override bool Equals(IIndexable2D other, double tolerance = 1e-13)
		{
			if ((this.NumRows != other.NumRows) || (this.NumColumns != other.NumColumns)) return false;
			var comparer = new ValueComparer(tolerance);
			for (int j = 0; j < NumColumns; ++j)
			{
				int colStart = colOffsets[j]; // Inclusive
				int colEnd = colOffsets[j + 1]; // Exclusive
				int previousRow = 0;
				for (int k = colStart; k < colEnd; ++k)
				{
					int row = rowIndices[k];
					for (int i = previousRow; i < row; ++i) // Zero entries between the stored ones
					{
						if (!comparer.AreEqual(0.0, other[i, j])) return false;
					}
					if (!comparer.AreEqual(values[k], other[row, j])) return false; // Non zero entry
					previousRow = row + 1;
				}
			}
			return true; // At this point all entries have been checked and are equal
		}

		public override Vector GetColumn(int colIndex)
		{
			Preconditions.CheckIndexCol(this, colIndex);
			double[] colVector = new double[NumRows];
			for (int k = colOffsets[colIndex]; k < colOffsets[colIndex + 1]; ++k) colVector[rowIndices[k]] = values[k];
			return Vector.CreateFromArray(colVector, false);
		}

		public override Vector GetRow(int rowIndex)
		{
			Preconditions.CheckIndexRow(this, rowIndex);
			double[] rowVector = new double[NumColumns];
			for (int j = 0; j < NumColumns; ++j)
			{
				int entryOffset = FindOffsetOf(rowIndex, j);
				if (entryOffset != zeroEntryOffset) rowVector[j] = values[entryOffset];
			}
			return Vector.CreateFromArray(rowVector, false);
		}

		public SparseFormat GetSparseFormat()
		{
			var format = new SparseFormat();
			format.RawValuesTitle = "Values";
			format.RawValuesArray = values;
			format.RawIndexArrays.Add("Row indices", rowIndices);
			format.RawIndexArrays.Add("Column offsets", colOffsets);
			return format;
		}

		public override bool HasSameFormat(CscMatrix other)
		{
			return (this.rowIndices == other.rowIndices) && (this.colOffsets == other.colOffsets);
		}

		public override Matrix MultiplyLeft(IReadOnlyMatrix other, bool transposeThis = false, bool transposeOther = false)
		{
			//TODO: To use BLAS for this too, we must accept row major matrices as output.
			if (transposeOther)
			{
				if (transposeThis)
				{
					Preconditions.CheckMultiplicationDimensions(other.NumRows, this.NumColumns);
					var result = Matrix.CreateZero(other.NumColumns, this.NumRows);
					CsrMultiplications.MatrixTransTimesCsr(this.NumColumns, values, colOffsets, rowIndices, other, result);
					return result;
				}
				else
				{
					Preconditions.CheckMultiplicationDimensions(other.NumRows, this.NumRows);
					var result = Matrix.CreateZero(other.NumColumns, this.NumColumns);
					CsrMultiplications.MatrixTransTimesCsrTrans(this.NumColumns, values, colOffsets, rowIndices, other, result);
					return result;
				}
			}
			else
			{
				if (transposeThis)
				{
					Preconditions.CheckMultiplicationDimensions(other.NumColumns, this.NumColumns);
					var result = Matrix.CreateZero(other.NumRows, this.NumRows);
					CsrMultiplications.MatrixTimesCsr(this.NumColumns, values, colOffsets, rowIndices, other, result);
					return result;
				}
				else
				{
					Preconditions.CheckMultiplicationDimensions(other.NumColumns, this.NumRows);
					var result = Matrix.CreateZero(other.NumRows, this.NumColumns);
					CsrMultiplications.MatrixTimesCsrTrans(this.NumColumns, values, colOffsets, rowIndices, other, result);
					return result;
				}
			}
		}

		public override Matrix MultiplyRight(IReadOnlyMatrix other, bool transposeThis = false, bool transposeOther = false)
		{
			// TODO: Throwing exceptions when csc is on the left seems attractive.
			if (transposeOther)
			{
				if (transposeThis)
				{
					Preconditions.CheckMultiplicationDimensions(this.NumRows, other.NumColumns);
					var result = Matrix.CreateZero(this.NumColumns, other.NumRows);
					CsrMultiplications.CsrTimesMatrixTrans(this.NumColumns, values, colOffsets, rowIndices, other, result);
					return result;
				}
				else
				{
					Preconditions.CheckMultiplicationDimensions(this.NumColumns, other.NumColumns);
					var result = Matrix.CreateZero(this.NumRows, other.NumRows);
					CsrMultiplications.CsrTransTimesMatrixTrans(this.NumColumns, values, colOffsets, rowIndices, other, result);
					return result;
				}
			}
			else
			{
				//TODO: perhaps I can use the left multiplications if the other matrix is also transposed
				if (other is Matrix dense) return MultiplyRight(dense, transposeThis);

				if (transposeThis)
				{
					Preconditions.CheckMultiplicationDimensions(this.NumRows, other.NumRows);
					var result = Matrix.CreateZero(this.NumColumns, other.NumColumns);
					CsrMultiplications.CsrTimesMatrix(this.NumColumns, values, colOffsets, rowIndices, other, result);
					return result;
				}
				else
				{
					Preconditions.CheckMultiplicationDimensions(this.NumColumns, other.NumRows);
					var result = Matrix.CreateZero(this.NumRows, other.NumColumns);
					CsrMultiplications.CsrTransTimesMatrix(this.NumColumns, values, colOffsets, rowIndices, other, result);
					return result;
				}
			}
		}

		/// <summary>
		/// Performs the matrix-matrix multiplication: oper(this) * <paramref name="other"/>.
		/// </summary>
		/// <param name="other">
		/// A matrix such that the <see cref="IIndexable2D.NumRows"/> of <paramref name="other"/> are equal to the 
		/// <see cref="IIndexable2D.NumColumns"/> of oper(this).
		/// </param>
		/// <param name="transposeThis">If true, oper(this) = transpose(this). Otherwise oper(this) = this.</param>
		/// <exception cref="NonMatchingDimensionsException">
		/// Thrown if <paramref name="otherMatrix"/> has different <see cref="IIndexable2D.NumRows"/> than the 
		/// <see cref="IIndexable2D.NumColumns"/> of oper(this).
		/// </exception>
		public Matrix MultiplyRight(Matrix other, bool transposeThis)
		{
			int numRowsResult;
			if (transposeThis)
			{
				Preconditions.CheckMultiplicationDimensions(this.NumRows, other.NumRows);
				numRowsResult = this.NumColumns;
			}
			else
			{
				Preconditions.CheckMultiplicationDimensions(this.NumColumns, other.NumRows);
				numRowsResult = this.NumRows;
			}

			var result = Matrix.CreateZero(numRowsResult, other.NumColumns);
			GlobalProvider.SparseBlas.Dcscgemm(transposeThis, this.NumRows, other.NumColumns, this.NumColumns, values, 
				colOffsets, rowIndices, other.RawData, result.RawData);
			return result;
		}

		public override IVector Multiply(IReadOnlyVector vector, bool transposeThis = false)
		{
			if (vector is Vector dense) return Multiply(dense, transposeThis);

			if (transposeThis)
			{
				var result = new double[NumColumns];
				Preconditions.CheckMultiplicationDimensions(NumRows, vector.Length);
				CsrMultiplications.CsrTimesVector(NumColumns, values, colOffsets, rowIndices, vector, result);
				return Vector.CreateFromArray(result, false);
			}
			else
			{
				var result = new double[NumRows];
				Preconditions.CheckMultiplicationDimensions(NumColumns, vector.Length);
				CsrMultiplications.CsrTransTimesVector(NumColumns, values, colOffsets, rowIndices, vector, result);
				return Vector.CreateFromArray(result, false);
			}
		}

		/// <summary>
		/// Performs the matrix-vector multiplication: oper(this) * <paramref name="vector"/>.
		/// To multiply this * columnVector, set <paramref name="transposeThis"/> to false.
		/// To multiply rowVector * this, set <paramref name="transposeThis"/> to true.
		/// </summary>
		/// <param name="vector">A vector with <see cref="IIndexable1D.Length"/> being equal to the 
		///     <see cref="IIndexable2D.NumColumns"/> of oper(this).</param>
		/// <param name="transposeThis">If true, oper(this) = transpose(this). Otherwise oper(this) = this.</param>
		/// <exception cref="NonMatchingDimensionsException">Thrown if the <see cref="IIndexable1D.Length"/> of
		///     <paramref name="vector"/> is different than the <see cref="NumColumns"/> of oper(this).</exception>
		public Vector Multiply(Vector vector, bool transposeThis = false)
		{
			//TODO: this performs redundant dimension checks, including checking the transposeThis flag.
			var result = Vector.CreateZero(transposeThis ? NumColumns : NumRows);
			MultiplyIntoResult(vector, result, transposeThis);
			return result;
		}

		public override void MultiplyIntoResult(IReadOnlyVector lhsVector, IVector rhsVector, bool transposeThis = false)
		{
			if (this.values.Length == 0)
			{
				Preconditions.CheckMultiplicationDimensions(NumColumns, lhsVector.Length);
				Preconditions.CheckSystemSolutionDimensions(NumRows, rhsVector.Length);
				return;
			}

			if ((lhsVector is Vector lhsDense) && (rhsVector is Vector rhsDense))
			{
				MultiplyIntoResult(lhsDense, rhsDense, transposeThis);
			}

			if (transposeThis)
			{
				Preconditions.CheckMultiplicationDimensions(NumRows, lhsVector.Length);
				Preconditions.CheckSystemSolutionDimensions(NumColumns, rhsVector.Length);
				CsrMultiplications.CsrTimesVector(NumColumns, values, colOffsets, rowIndices, lhsVector, rhsVector);
			}
			else
			{
				Preconditions.CheckMultiplicationDimensions(NumColumns, lhsVector.Length);
				Preconditions.CheckSystemSolutionDimensions(NumRows, rhsVector.Length);
				CsrMultiplications.CsrTransTimesVector(NumColumns, values, colOffsets, rowIndices, lhsVector, rhsVector);
			}
		}

		/// <summary>
		/// Performs the matrix-vector multiplication: <paramref name="rhsVector"/> = oper(this) * <paramref name="vector"/>.
		/// To multiply this * columnVector, set <paramref name="transposeThis"/> to false.
		/// To multiply rowVector * this, set <paramref name="transposeThis"/> to true.
		/// The resulting vector will overwrite the entries of <paramref name="rhsVector"/>.
		/// </summary>
		/// <param name="lhsVector">
		/// The vector that will be multiplied by this matrix. It sits on the left hand side of the equation y = oper(A) * x.
		/// Constraints: <paramref name="lhsVector"/>.<see cref="IIndexable1D.Length"/> 
		/// == oper(this).<see cref="IIndexable2D.NumColumns"/>.
		/// </param>
		/// <param name="rhsVector">
		/// The vector that will be overwritten by the result of the multiplication. It sits on the right hand side of the 
		/// equation y = oper(A) * x. Constraints: <paramref name="lhsVector"/>.<see cref="IIndexable1D.Length"/> 
		/// == oper(this).<see cref="IIndexable2D.NumRows"/>.
		/// </param>
		/// <param name="transposeThis">If true, oper(this) = transpose(this). Otherwise oper(this) = this.</param>
		/// <exception cref="NonMatchingDimensionsException">
		/// Thrown if the <see cref="IIndexable1D.Length"/> of <paramref name="lhsVector"/> or <paramref name="rhsVector"/> 
		/// violate the described contraints.
		/// </exception>
		public void MultiplyIntoResult(Vector lhsVector, Vector rhsVector, bool transposeThis = false)
		{
			if (this.values.Length == 0)
			{
				Preconditions.CheckMultiplicationDimensions(NumColumns, lhsVector.Length);
				Preconditions.CheckSystemSolutionDimensions(NumRows, rhsVector.Length);
				return;
			}

			if (transposeThis)
			{
				Preconditions.CheckMultiplicationDimensions(NumRows, lhsVector.Length);
				Preconditions.CheckSystemSolutionDimensions(NumColumns, rhsVector.Length);
			}
			else
			{
				Preconditions.CheckMultiplicationDimensions(NumColumns, lhsVector.Length);
				Preconditions.CheckSystemSolutionDimensions(NumRows, rhsVector.Length);
			}

			GlobalProvider.SparseBlas.Dcscgemv(transposeThis, NumRows, NumColumns, values, colOffsets, rowIndices,
					lhsVector.RawData, 0, rhsVector.RawData, 0);
		}

		public override double Reduce(
			double identityValue, ProcessEntry processEntry, ProcessZeros processZeros, Finalize finalize)
			=> ReduceNonSymmetric(identityValue, processEntry, processZeros, finalize);

		public override IMatrix Transpose() => TransposeToCSR(true);

		/// <summary>
		/// Creates a new <see cref="CscMatrix"/> instance, that is transpose to this: result[i, j] = this[j, i].
		/// </summary>
		public CscMatrix TransposeToCSC()
		{
			// Use C# port of the scipy method.
			// TODO: Perhaps it could be done faster by making extra assumptions. Otherwise use SparseBLAS
			int nnz = this.values.Length;
			var csrValues = new double[nnz];
			var csrColIndices = new int[nnz];
			var csrRowOffsets = new int[NumRows + 1];

			Conversions.CsrToCsc(NumColumns, NumRows, this.colOffsets, this.rowIndices, this.values,
				csrRowOffsets, csrColIndices, csrValues);

			return new CscMatrix(NumColumns, NumRows, csrValues, csrColIndices, csrRowOffsets);
		}

		/// <summary>
		/// Creates a new <see cref="CsrMatrix"/> instance, that is transpose to this: result[i, j] = this[j, i]. The 
		/// internal arrays can be copied or shared with this <see cref="CscMatrix"/> instance.
		/// </summary>
		/// <param name="copyInternalArray">If true, the internal arrays that store the entries of this 
		///     <see cref="CscMatrix"/> instance will be copied and the new <see cref="CsrMatrix"/> instance 
		///     instance will have references to the copies, which is safer. If false, both the new matrix and this one will have  
		///     references to the same internal arrays, which is faster.</param>
		public CsrMatrix TransposeToCSR(bool copyInternalArrays)
		{
			if (copyInternalArrays)
			{
				double[] valuesCopy = new double[values.Length];
				Array.Copy(values, valuesCopy, values.Length);
				int[] rowIndicesCopy = new int[rowIndices.Length];
				Array.Copy(rowIndices, rowIndicesCopy, rowIndices.Length);
				int[] colOffsetsCopy = new int[colOffsets.Length];
				Array.Copy(colOffsets, colOffsetsCopy, colOffsets.Length);
				return CsrMatrix.CreateFromArrays(NumColumns, NumRows, valuesCopy, rowIndicesCopy, colOffsetsCopy, false);
			}
			else return CsrMatrix.CreateFromArrays(NumColumns, NumRows, values, rowIndices, colOffsets, false);
		}

		/// <summary>
		/// Return the index into values and rowIndices arrays, if the (rowIdx, colIdx) entry is within the pattern. 
		/// Otherwise returns <see cref="zeroEntryOffset"/>.
		/// </summary>
		/// <param name="rowIdx"></param>
		/// <param name="colIdx"></param>
		private int FindOffsetOf(int rowIdx, int colIdx)
		{
			//TODO: if we have a true bool flag to indicate that the row indices of each column are sorted, then use binary search.
			Preconditions.CheckIndices(this, rowIdx, colIdx); //TODO: check indices?
			int colStart = colOffsets[colIdx]; //inclusive
			int colEnd = colOffsets[colIdx + 1]; //exclusive
			for (int k = colStart; k < colEnd; ++k)
			{
				if (rowIndices[k] == rowIdx) return k;
			}
			return zeroEntryOffset;
		}
	}
}
