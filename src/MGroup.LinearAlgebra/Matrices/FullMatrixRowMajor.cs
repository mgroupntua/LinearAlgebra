namespace MGroup.LinearAlgebra.Matrices
{
	using System;
	using System.Collections.Generic;
	using System.Diagnostics;
	using System.Text;

	using MGroup.LinearAlgebra.Commons;
	using MGroup.LinearAlgebra.Implementations;
	using MGroup.LinearAlgebra.Reduction;
	using MGroup.LinearAlgebra.Vectors;

	using static MGroup.LinearAlgebra.LibrarySettings;

	public class FullMatrixRowMajor : IMatrixView
	{
		private readonly double[] values;

		private FullMatrixRowMajor(int numRows, int numColumns, double[] values)
		{
			NumRows = numRows;
			NumColumns = numColumns;
			this.values = values;
		}

		public double this[int rowIdx, int colIdx]
		{
			get => values[rowIdx * NumColumns + colIdx];
			set => values[rowIdx * NumColumns + colIdx] = value;
		}

		public int NumColumns { get; }

		public int NumRows { get; }

		public double[] RawData => values;

		public MatrixSymmetry MatrixSymmetry => MatrixSymmetry.Unknown;

		public static FullMatrixRowMajor CreateFromArray(int numRows, int numColumns, double[] values, bool copyArray = false)
		{
			if (copyArray)
			{
				var clone = new double[values.Length];
				Array.Copy(values, clone, clone.Length);
				return new FullMatrixRowMajor(numRows, numColumns, clone);
			}
			else
			{
				return new FullMatrixRowMajor(numRows, numColumns, values);
			}
		}

		/// <summary>
		/// Initializes a new instance of <see cref="FullMatrixRowMajor"/> by copying the entries of <paramref name="array2D"/>.
		/// </summary>
		/// <param name="array2D">A 2-dimensional array containing the entries of the matrix. It will be copied.</param>
		public static FullMatrixRowMajor CreateFromArray(double[,] array2D)
		{
			int numRows = array2D.GetLength(0);
			int numCols = array2D.GetLength(1);
			return new FullMatrixRowMajor(numRows, numCols, Conversions.Array2DToFullRowMajor(array2D));
		}

		public static FullMatrixRowMajor CreateZero(int numRows, int numColumns)
		{
			return new FullMatrixRowMajor(numRows, numColumns, new double[numRows * numColumns]);
		}

		public IMatrix Axpy(IMatrixView otherMatrix, double otherCoefficient)
		{
			throw new NotImplementedException();
		}

		public IMatrix Copy(bool copyIndexingData = false)
		{
			throw new NotImplementedException();
		}

		public Matrix CopyToFullMatrix()
		{
			throw new NotImplementedException();
		}

		public IMatrix DoEntrywise(IMatrixView matrix, Func<double, double, double> binaryOperation)
		{
			throw new NotImplementedException();
		}

		public IMatrix DoToAllEntries(Func<double, double> unaryOperation)
		{
			throw new NotImplementedException();
		}

		public bool Equals(IIndexable2D other, double tolerance = 1E-13)
		{
			throw new NotImplementedException();
		}

		public Vector GetColumn(int colIndex)
		{
			throw new NotImplementedException();
		}

		public Vector GetRow(int rowIndex)
		{
			throw new NotImplementedException();
		}

		public IMatrix GetSubmatrix(int[] rowIndices, int[] colIndices)
		{
			throw new NotImplementedException();
		}

		public IMatrix GetSubmatrix(int rowStartInclusive, int rowEndExclusive, int colStartInclusive, int colEndExclusive)
		{
			throw new NotImplementedException();
		}

		public IMatrix LinearCombination(double thisCoefficient, IMatrixView otherMatrix, double otherCoefficient)
		{
			throw new NotImplementedException();
		}

		public IVector Multiply(IVectorView vector, bool transposeThis = false)
		{
			if (vector is Vector dense) return Multiply(dense, transposeThis);
			else throw new NotImplementedException();
		}

		public Vector Multiply(Vector vector, bool transposeThis = false)
		{
			//TODO: this performs redundant dimension checks, including checking the transposeThis flag.
			var result = Vector.CreateZero(transposeThis ? NumColumns : NumRows);
			MultiplyIntoResult(vector, result, transposeThis);
			return result;
		}

		public void MultiplyIntoResult(IVectorView lhsVector, IVector rhsVector, bool transposeThis = false)
		{
			if (lhsVector is Vector denseLhs && rhsVector is Vector denseRhs)
			{
				MultiplyIntoResult(denseLhs, denseRhs, transposeThis);
			}
			else
			{
				throw new NotImplementedException();
			}
		}

		public void MultiplyIntoResult(Vector lhsVector, Vector rhsVector, bool transposeThis = false)
		{
			(TransposeMatrix transposeA, int lhsLength, int rhsLength) = TransposeUtilities.PrepareBlas(this, transposeThis);
			Preconditions.CheckMultiplicationDimensions(lhsLength, lhsVector.Length);
			Preconditions.CheckSystemSolutionDimensions(rhsLength, rhsVector.Length);
			Blas.DgemvRowMajor(transposeA, NumRows, NumColumns, this.values, lhsVector.RawData, rhsVector.RawData);
		}

		public Matrix MultiplyLeft(IMatrixView other, bool transposeThis = false, bool transposeOther = false)
		{
			throw new NotImplementedException();
		}

		public Matrix MultiplyRight(IMatrixView other, bool transposeThis = false, bool transposeOther = false)
		{
			throw new NotImplementedException();
		}

		public double Reduce(double identityValue, ProcessEntry processEntry, ProcessZeros processZeros, Finalize finalize)
		{
			throw new NotImplementedException();
		}

		public void SetRow(int rowIdx, Vector rowValues)
		{
			Debug.Assert(rowIdx >= 0 && rowIdx < NumRows);
			Debug.Assert(rowValues.Length == NumColumns);
			Array.Copy(rowValues.RawData, 0, values, rowIdx * NumColumns, NumColumns);
		}

		public IMatrix Scale(double scalar)
		{
			throw new NotImplementedException();
		}

		public IMatrix Transpose()
		{
			throw new NotImplementedException();
		}
	}
}
