namespace MGroup.LinearAlgebra.Matrices
{
	using System;
	using System.Diagnostics;

	using MGroup.LinearAlgebra.Commons;
	using MGroup.LinearAlgebra.Implementations;
	using MGroup.LinearAlgebra.Reduction;
	using MGroup.LinearAlgebra.Vectors;

	using static MGroup.LinearAlgebra.LibrarySettings;

	public sealed class FullMatrixRowMajor : DefaultMatrix
	{
		private readonly double[] values;

		private FullMatrixRowMajor(int numRows, int numColumns, double[] values)
		{
			NumRows = numRows;
			NumColumns = numColumns;
			this.values = values;
		}

		public override int NumColumns { get; }

		public override int NumRows { get; }

		public double[] RawData => values;

		public override double this[int rowIdx, int colIdx]
		{
			get => values[rowIdx * NumColumns + colIdx];
			set => values[rowIdx * NumColumns + colIdx] = value;
		}

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

		public override void Clear() => Array.Clear(values, 0, values.Length);

		public override IMatrix CreateZeroMatrixWithSameFormat()
			=> new FullMatrixRowMajor(NumRows, NumColumns, new double[values.Length]);

		public override bool HasSameFormat(IMatrixView other)
		{
			if (other is FullMatrixRowMajor casted)
			{
				return (this.NumRows == casted.NumRows) && (this.NumColumns == other.NumColumns);
			}

			return false;
		}

		public override IVector Multiply(IVectorView vector, bool transposeThis = false)
		{
			if (vector is Vector dense)
			{
				return Multiply(dense, transposeThis);
			}
			else
			{
				return base.Multiply(vector, transposeThis);
			}
		}

		public Vector Multiply(Vector vector, bool transposeThis = false)
		{
			//TODO: this performs redundant dimension checks, including checking the transposeThis flag.
			var result = Vector.CreateZero(transposeThis ? NumColumns : NumRows);
			MultiplyIntoResult(vector, result, transposeThis);
			return result;
		}

		public override void MultiplyIntoResult(IVectorView lhsVector, IVector rhsVector, bool transposeThis = false)
		{
			if (lhsVector is Vector denseLhs && rhsVector is Vector denseRhs)
			{
				MultiplyIntoResult(denseLhs, denseRhs, transposeThis);
			}
			else
			{
				base.MultiplyIntoResult(lhsVector, rhsVector, transposeThis);
			}
		}

		public void MultiplyIntoResult(Vector lhsVector, Vector rhsVector, bool transposeThis = false)
		{
			(TransposeMatrix transposeA, int lhsLength, int rhsLength) = TransposeUtilities.PrepareBlas(this, transposeThis);
			Preconditions.CheckMultiplicationDimensions(lhsLength, lhsVector.Length);
			Preconditions.CheckSystemSolutionDimensions(rhsLength, rhsVector.Length);
			GlobalProvider.Blas.DgemvRowMajor(transposeA, NumRows, NumColumns, this.values, lhsVector.RawData, rhsVector.RawData);
		}

		public void SetRow(int rowIdx, Vector rowValues)
		{
			Debug.Assert(rowIdx >= 0 && rowIdx < NumRows);
			Debug.Assert(rowValues.Length == NumColumns);
			Array.Copy(rowValues.RawData, 0, values, rowIdx * NumColumns, NumColumns);
		}

		public override IMatrix Transpose() => Matrix.CreateFromArray(values, NumColumns, NumRows, true);
	}
}
