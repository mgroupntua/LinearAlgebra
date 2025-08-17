namespace MGroup.LinearAlgebra.Matrices
{
	using System;
	using System.Diagnostics;

	using MGroup.LinearAlgebra.Commons;
	using MGroup.LinearAlgebra.Exceptions;
	using MGroup.LinearAlgebra.Implementations;
	using MGroup.LinearAlgebra.Reduction;
	using MGroup.LinearAlgebra.Vectors;

	using static MGroup.LinearAlgebra.Commons.PerformanceWarnings;

	[Serializable]
	public abstract class DefaultMatrix : IMatrix
	{
		public virtual MatrixSymmetry MatrixSymmetry { get; } = MatrixSymmetry.Unknown;

		public abstract int NumColumns { get; }

		public abstract int NumRows { get; }

		/// <summary>
		/// The entry with row index = rowIdx and column index = colIdx.
		/// </summary>
		/// <param name="rowIdx">The row index: 0 &lt;= <paramref name="rowIdx"/> &lt; <see cref="NumRows"/>.</param>
		/// <param name="colIdx">The column index: 0 &lt;= <paramref name="colIdx"/> &lt; <see cref="NumColumns"/>.</param>
		/// <exception cref="IndexOutOfRangeException">
		/// Thrown if <paramref name="rowIdx"/> or <paramref name="colIdx"/> violate the described constraints.
		/// </exception>
		/// <exception cref="SparsityPatternModifiedException">
		/// Thrown if you try to set an entry that is not explicitly stored by the matrix storage format
		/// (e.g. structural zeros in sparse formats).
		/// </exception>
		/// <remarks>
		/// Indexing is inefficient for most matrix storage formats. There is usually a method to do any job more efficiently.
		/// </remarks>
		public abstract double this[int rowIdx, int colIdx] { get; set; }

		public abstract void Clear();

		public abstract IMatrix CreateZeroMatrixWithSameFormat();

		public abstract bool HasSameFormat(IMatrixView otherMatrix);

		public virtual IMatrix Axpy(IMatrixView otherMatrix, double otherCoefficient)
			=> LinearCombination(1.0, otherMatrix, otherCoefficient);

		public virtual void AxpyIntoThis(IMatrixView otherMatrix, double otherCoefficient)
			=> LinearCombinationIntoThis(1.0, otherMatrix, otherCoefficient);

		public virtual IMatrix Copy(bool copyIndexingData = false)
		{
			WarnAboutPerformanceBottlenecks();
			ProhibitPerformanceBottlenecks();
			IMatrix clone = CreateZeroMatrixWithSameFormat();
			for (int j = 0; j < NumColumns; ++j)
			{
				for (int i = 0; i < NumRows; ++i)
				{
					double val = this[i, j];
					if (val != 0.0)
					{
						clone.Set(i, j, val);
					}
				}
			}

			return clone;
		}

		public virtual Matrix CopyToFullMatrix()
		{
			WarnAboutPerformanceBottlenecks();
			ProhibitPerformanceBottlenecks();
			var result = Matrix.CreateZero(NumRows, NumColumns);
			for (int j = 0; j < NumColumns; ++j)
			{
				for (int i = 0; i < NumRows; ++i)
				{
					result[i, j] = this[i, j];
				}
			}

			return result;
		}

		public virtual IMatrix DoEntrywise(IMatrixView otherMatrix, Func<double, double, double> binaryOperation)
		{
			Preconditions.CheckSameMatrixDimensions(this, otherMatrix);
			if (this.HasSameFormat(otherMatrix) && (binaryOperation(0.0, 0.0) == 0.0))
			{
				IMatrix result = Copy(copyIndexingData: false);
				result.DoEntrywiseIntoThis(otherMatrix, binaryOperation);
				return result;
			}
			else
			{
				WarnAboutPerformanceBottlenecks();
				ProhibitPerformanceBottlenecks();
				var result = Matrix.CreateZero(NumRows, NumColumns);
				for (int j = 0; j < NumColumns; ++j)
				{
					for (int i = 0; i < NumRows; ++i)
					{
						result[i, j] = binaryOperation(this[i, j], otherMatrix[i, j]);
					}
				}

				return result;
			}
		}

		public virtual void DoEntrywiseIntoThis(IMatrixView otherMatrix, Func<double, double, double> binaryOperation)
		{
			Preconditions.CheckSameMatrixDimensions(this, otherMatrix);
			WarnAboutPerformanceBottlenecks();
			ProhibitPerformanceBottlenecks();
			for (int j = 0; j < NumColumns; ++j)
			{
				for (int i = 0; i < NumRows; ++i)
				{
					this[i, j] = binaryOperation(this[i, j], otherMatrix[i, j]);
				}
			}
		}

		public virtual IMatrix DoToAllEntries(Func<double, double> unaryOperation)
		{
			WarnAboutPerformanceBottlenecks();
			ProhibitPerformanceBottlenecks();
			IMatrix result;
			if (unaryOperation(0.0) == 0.0)
			{
				result = Copy(copyIndexingData: false);
			}
			else
			{
				result = CopyToFullMatrix();
			}

			result.DoToAllEntriesIntoThis(unaryOperation);
			return result;
		}

		public virtual void DoToAllEntriesIntoThis(Func<double, double> unaryOperation)
		{
			WarnAboutPerformanceBottlenecks();
			ProhibitPerformanceBottlenecks();
			for (int j = 0; j < NumColumns; ++j)
			{
				for (int i = 0; i < NumRows; ++i)
				{
					this[i, j] = unaryOperation(this[i, j]);
				}
			}
		}

		public virtual bool Equals(IIndexable2D other, double tolerance = 1E-13)
		{
			if ((this.NumRows != other.NumRows) || (this.NumColumns != other.NumColumns))
			{
				return false;
			}

			WarnAboutPerformanceBottlenecks();
			ProhibitPerformanceBottlenecks();

			var comparer = new ValueComparer(tolerance);
			for (int j = 0; j < NumColumns; ++j)
			{
				for (int i = 0; i < NumRows; ++i)
				{
					if (!comparer.AreEqual(this[i, j], other[i, j]))
					{
						return false;
					}
				}
			}

			return true;
		}

		public virtual Vector GetDiagonal() => Vector.CreateFromArray(this.GetDiagonalAsArray(), false);

		public virtual double[] GetDiagonalAsArray()
		{
			Preconditions.CheckSquare(this);
			WarnAboutPerformanceBottlenecks();
			ProhibitPerformanceBottlenecks();
			double[] diag = new double[NumRows];
			for (int i = 0; i < NumRows; ++i)
			{
				diag[i] = this[i, i];
			}

			return diag;
		}

		public virtual Vector GetColumn(int colIndex)
		{
			Preconditions.CheckIndexCol(this, colIndex);
			WarnAboutPerformanceBottlenecks();
			ProhibitPerformanceBottlenecks();
			var columnVector = new double[NumRows];
			for (int i = 0; i < NumRows; ++i)
			{
				columnVector[i] = this[i, colIndex];
			}

			return Vector.CreateFromArray(columnVector, false);
		}

		public virtual Vector GetRow(int rowIndex)
		{
			WarnAboutPerformanceBottlenecks();
			ProhibitPerformanceBottlenecks();
			Preconditions.CheckIndexRow(this, rowIndex);
			var rowVector = new double[NumColumns];
			for (int j = 0; j < NumColumns; ++j)
			{
				rowVector[j] = this[rowIndex, j];
			}

			return Vector.CreateFromArray(rowVector, false);
		}

		public virtual IMatrix GetSubmatrix(int[] rowIndices, int[] colIndices)
		{
			WarnAboutPerformanceBottlenecks();
			ProhibitPerformanceBottlenecks();
			var submatrix = Matrix.CreateZero(rowIndices.Length, colIndices.Length);
			for (int j = 0; j < colIndices.Length; ++j)
			{
				for (int i = 0; i < rowIndices.Length; ++i)
				{
					submatrix[i, j] = this[rowIndices[i], colIndices[j]];
				}
			}

			return submatrix;
		}

		public virtual IMatrix GetSubmatrix(int rowStartInclusive, int rowEndExclusive, int colStartInclusive, int colEndExclusive)
		{
			if (rowStartInclusive < 0 || rowEndExclusive >= NumRows || colStartInclusive < 0 || colEndExclusive >= NumColumns)
			{
				throw new NonMatchingDimensionsException(
					"The submatrix cannot extend outside the bounds of the original matrix.");
			}

			WarnAboutPerformanceBottlenecks();
			ProhibitPerformanceBottlenecks();
			int numNewRows = rowEndExclusive - rowStartInclusive;
			int numNewCols = colEndExclusive - colStartInclusive;
			var submatrix = Matrix.CreateZero(numNewRows, numNewCols);
			for (int j = 0; j < numNewCols; ++j)
			{
				for (int i = 0; i < numNewRows; ++i)
				{
					submatrix[i, j] = this[rowStartInclusive + i, colStartInclusive + j];
				}
			}

			return submatrix;
		}

		public virtual IMatrix LinearCombination(double thisCoefficient, IMatrixView otherMatrix, double otherCoefficient)
			=> DoEntrywise(otherMatrix, (x, y) => thisCoefficient * x + otherCoefficient * y);

		public virtual void LinearCombinationIntoThis(double thisCoefficient, IMatrixView otherMatrix, double otherCoefficient)
			=> DoEntrywiseIntoThis(otherMatrix, (x, y) => thisCoefficient * x + otherCoefficient * y);

		public virtual IVector Multiply(IVectorView vector, bool transposeThis = false)
		{
			Vector result = transposeThis ? Vector.CreateZero(this.NumColumns) : Vector.CreateZero(this.NumRows);
			MultiplyIntoResult(vector, result, transposeThis);
			return result;
		}

		public virtual void MultiplyIntoResult(IVectorView lhsVector, IVector rhsVector, bool transposeThis = false)
		{
			WarnAboutPerformanceBottlenecks();
			ProhibitPerformanceBottlenecks();
			if (transposeThis)
			{
				Preconditions.CheckMultiplicationDimensionsMatrixVectorTranspose(this, lhsVector, rhsVector);
				for (int i = 0; i < rhsVector.Length; ++i)
				{
					for (int j = 0; j < lhsVector.Length; ++j)
					{
						rhsVector.Set(i, this[j, i] * lhsVector[j]);
					}
				}
			}
			else
			{
				Preconditions.CheckMultiplicationDimensionsMatrixVector(this, lhsVector, rhsVector);
				for (int i = 0; i < rhsVector.Length; ++i)
				{
					for (int j = 0; j < lhsVector.Length; ++j)
					{
						rhsVector.Set(i, this[i, j] * lhsVector[j]);
					}
				}
			}
		}

		public virtual Matrix MultiplyLeft(IMatrixView otherMatrix, bool transposeThis = false, bool transposeOther = false)
			=> MultiplyMatrices(otherMatrix, this, transposeOther, transposeThis);

		public virtual Matrix MultiplyRight(IMatrixView otherMatrix, bool transposeThis = false, bool transposeOther = false)
			=> MultiplyMatrices(this, otherMatrix, transposeThis, transposeOther);

		public virtual double Reduce(
			double identityValue, ProcessEntry processEntry, ProcessZeros processZeros, Finalize finalize)
		{
			WarnAboutPerformanceBottlenecks();
			ProhibitPerformanceBottlenecks();
			double accumulator = identityValue;
			for (int j = 0; j < NumColumns; ++j)
			{
				for (int i = 0; i < NumRows; ++i)
				{
					accumulator = processEntry(this[i, j], accumulator);
				}
			}

			return finalize(accumulator);
		}

		public virtual IMatrix Scale(double scalar) => DoToAllEntries(x => scalar * x);

		public virtual void ScaleIntoThis(double scalar) => DoToAllEntriesIntoThis(x => scalar * x);

		public virtual void Set(int rowIdx, int colIdx, double value) => this[rowIdx, colIdx] = value;

		public virtual IMatrix Transpose()
		{
			if (MatrixSymmetry == MatrixSymmetry.Symmetric)
			{
				return Copy();
			}

			WarnAboutPerformanceBottlenecks();
			ProhibitPerformanceBottlenecks();
			var result = Matrix.CreateZero(NumRows, NumColumns);
			for (int i = 0; i < NumRows; ++i)
			{
				for (int j = 0; j < NumColumns; ++j)
				{
					result[j, i] = this[i, j];
				}
			}

			return result;
		}

		private static Matrix MultiplyMatrices(
			IMatrixView matrixLeft, IMatrixView matrixRight, bool transposeLeft, bool transposeRight)
		{
			WarnAboutPerformanceBottlenecks();
			ProhibitPerformanceBottlenecks();
			if (transposeLeft)
			{
				if (transposeRight)
				{
					Preconditions.CheckMultiplicationDimensions(matrixLeft.NumRows, matrixRight.NumColumns);
					var result = Matrix.CreateZero(matrixLeft.NumColumns, matrixRight.NumRows);
					for (int i = 0; i < result.NumRows; ++i)
					{
						for (int j = 0; j < result.NumColumns; ++j)
						{
							for (int k = 0; k < matrixLeft.NumRows; ++k)
							{
								result[i, j] += matrixLeft[k, i] * matrixRight[j, k];
							}
						}
					}

					return result;
				}
				else
				{
					Preconditions.CheckMultiplicationDimensions(matrixLeft.NumRows, matrixRight.NumRows);
					var result = Matrix.CreateZero(matrixLeft.NumColumns, matrixRight.NumColumns);
					for (int i = 0; i < result.NumRows; ++i)
					{
						for (int j = 0; j < result.NumColumns; ++j)
						{
							for (int k = 0; k < matrixLeft.NumRows; ++k)
							{
								result[i, j] += matrixLeft[k, i] * matrixRight[k, j];
							}
						}
					}

					return result;
				}
			}
			else
			{
				if (transposeRight)
				{
					Preconditions.CheckMultiplicationDimensions(matrixLeft.NumColumns, matrixRight.NumColumns);
					var result = Matrix.CreateZero(matrixLeft.NumRows, matrixRight.NumRows);
					for (int i = 0; i < result.NumRows; ++i)
					{
						for (int j = 0; j < result.NumColumns; ++j)
						{
							for (int k = 0; k < matrixLeft.NumColumns; ++k)
							{
								result[i, j] += matrixLeft[i, k] * matrixRight[j, k];
							}
						}
					}

					return result;
				}
				else
				{
					Preconditions.CheckMultiplicationDimensions(matrixLeft.NumColumns, matrixRight.NumRows);
					var result = Matrix.CreateZero(matrixLeft.NumRows, matrixRight.NumColumns);
					for (int i = 0; i < result.NumRows; ++i)
					{
						for (int j = 0; j < result.NumColumns; ++j)
						{
							for (int k = 0; k < matrixLeft.NumColumns; ++k)
							{
								result[i, j] += matrixLeft[i, k] * matrixRight[k, j];
							}
						}
					}

					return result;
				}
			}
		}
	}
}
