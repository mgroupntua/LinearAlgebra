//TODO: DoEntryWise may change the sparsity pattern if binaryOperation(0, 0)!=0. Fix this in all triangular and sparse matrices.
namespace MGroup.LinearAlgebra.Matrices
{
	using System;
	using System.Runtime.CompilerServices;

	using MGroup.LinearAlgebra.Commons;
	using MGroup.LinearAlgebra.Exceptions;
	using MGroup.LinearAlgebra.Implementations;
	using MGroup.LinearAlgebra.Reduction;
	using MGroup.LinearAlgebra.Vectors;

	using static MGroup.LinearAlgebra.LibrarySettings;

	/// <summary>
	/// Lower triangular square matrix in row major Packed storage format (only stores the n*(n+1)/2 non zeros). Uses LAPACK. 
	/// For the more information about the layout, see 
	/// <see cref="https://software.intel.com/en-us/mkl-developer-reference-c-matrix-storage-schemes-for-lapack-routines."/>
	/// </summary>
	public sealed class TriangularLower : ValuesBackedMatrix<TriangularLower>
	{
		/// <summary>
		/// Packed storage, row major order: L[i, j] = data[j + (i+1)*i/2] for 0 &lt;= j &lt;= i &lt; n.
		/// </summary>
		private readonly double[] values;

		private TriangularLower(double[] data, int order)
		{
			this.values = data;
			this.Order = order;
		}

		public override MatrixSymmetry MatrixSymmetry => MatrixSymmetry.NonSymmetric;

		public override int NumColumns => Order;

		public override int NumRows => Order;

		/// <summary>
		/// The number of rows/columns of the square matrix.
		/// </summary>
		public int Order { get; }

		/// <summary>
		/// The internal array that stores the entries of the lower triangle (packed storage format) in row major layout.
		/// It should only be used for passing the raw array to linear algebra libraries.
		/// </summary>
		public override double[] RawValues => values;

		public override double this[int rowIdx, int colIdx]
		{
			get
			{
				Preconditions.CheckIndices(this, rowIdx, colIdx);
				return (rowIdx >= colIdx) ? values[FindIndex1D(rowIdx, colIdx)] : 0.0;
			}

			set
			{
				Preconditions.CheckIndices(this, rowIdx, colIdx);
				if (rowIdx >= colIdx)
				{
					values[FindIndex1D(rowIdx, colIdx)] = value;
				}
				else
				{
					throw new SparsityPatternModifiedException(
						$"Cannot change the superdiagonal entry A[{rowIdx}, {colIdx}] = 0.");
				}
			}
		}

		/// <summary>
		/// Initializes a new instance of <see cref="TriangularLower"/> by copying the lower triangle of the provided 2D array.
		/// </summary>
		/// <param name="array2D">A 2-dimensional array containing the elements of the matrix. Constraints: 
		///     <paramref name="array2D"/>.GetLength(0) == <paramref name="array2D"/>.GetLength(1).</param>
		/// <exception cref="NonMatchingDimensionsException">Thrown if <paramref name="array2D"/>.GetLength(0) != 
		///     <paramref name="array2D"/>.GetLength(1).</exception>
		public static TriangularLower CreateFromArray(double[,] array2D)
		{
			int numRows = array2D.GetLength(0);
			int numCols = array2D.GetLength(1);
			if (numRows != numCols)
			{
				string msg = $"Provided array must have the same dimensions, but was ({numRows}x{numCols})";
				throw new NonMatchingDimensionsException(msg);
			}
			return new TriangularLower(Conversions.Array2DToPackedLowerRowMajor(array2D), numRows);
		}

		/// <summary>
		/// Initializes a new instance of <see cref="TriangularLower"/> with <paramref name="array1D"/> or a clone as its 
		/// internal array.
		/// </summary>
		/// <param name="order">The number of rows/columns of the new square matrix.</param>
		/// <param name="array1D">An 1-dimensional array containing the elements of the lower triangle of the matrix in row 
		///     major order.</param>
		/// <param name="copyArray">If true, <paramref name="array1D"/> will be copied and the new <see cref="TriangularLower"/>  
		///     instance will have a reference to the copy, which is safer. If false, the new matrix will have a reference to 
		///     <paramref name="array1D"/> itself, which is faster.</param>
		public static TriangularLower CreateFromArray(int order, double[] array1D, bool copyArray = false)
		{
			if (copyArray)
			{
				var clone = new double[array1D.Length];
				Array.Copy(array1D, clone, array1D.Length);
				return new TriangularLower(clone, order);
			}
			else return new TriangularLower(array1D, order);
		}

		/// <summary>
		/// Initializes a new instance of <see cref="TriangularLower"/> with all entries being equal to 0. Only the lower 
		/// triangle zero entries are explictily stored though.
		/// </summary>
		/// <param name="order">The number of rows/columns of the new square matrix.</param>
		public static TriangularLower CreateZero(int order)
		{
			var data = new double[(order * (order + 1)) / 2];
			return new TriangularLower(data, order);
		}

		/// <summary>
		/// Calculates the determinant of this matrix.
		/// </summary>
		public double CalcDeterminant()
		{
			// TODO: Find more effienct formulas for the diagonal accesses.
			double det = 1.0;
			for (int i = 0; i < Order; ++i)
			{
				det *= values[FindIndex1D(i, i)];
			}
			return det;
		}

		/// <summary>
		/// Copies the entries of this matrix.
		/// </summary>
		public override TriangularLower CopyAsSameType(bool copyIndexingData = false)
		{
			var clone = new double[this.values.Length];
			Array.Copy(this.values, clone, this.values.Length);
			return new TriangularLower(clone, Order);
		}

		/// <summary>
		/// Copies the entries of the matrix into a 2-dimensional array. The returned array has
		/// length(0) = length(1) = <see cref="Order"/>.
		/// </summary>
		public double[,] CopyToArray2D() => Conversions.PackedLowerRowMajorToArray2D(values);

		public override Matrix CopyToFullMatrix()
		{
			// TODO: This won't work if the implementation of Matrix changes
			double[] fullArray = Conversions.PackedLowerRowMajorToFullColMajor(values, Order);
			return Matrix.CreateFromArray(fullArray, Order, Order, false);
		}

		public override TriangularLower CreateZeroMatrixSame()
		{
			var resultValues = new double[values.Length];
			return new TriangularLower(resultValues, Order);
		}

		public override Vector GetColumn(int colIndex)
		{
			Preconditions.CheckIndexCol(this, colIndex);
			var columnVector = new double[NumRows];
			for (int i = colIndex; i < NumRows; ++i)
			{
				columnVector[i] = values[colIndex + (i * (i + 1)) / 2];
			}

			return Vector.CreateFromArray(columnVector);
		}

		public override Vector GetRow(int rowIndex)
		{
			Preconditions.CheckIndexRow(this, rowIndex);
			var rowVector = new double[NumColumns];
			int numNonZerosRow = rowIndex + 1;
			int rowOffset = (rowIndex * (rowIndex + 1)) / 2;
			Array.Copy(values, rowOffset, rowVector, 0, numNonZerosRow);
			return Vector.CreateFromArray(rowVector);
		}

		public override bool HasSameFormat(TriangularLower other) => this.Order == other.Order;

		public override IVector Multiply(IVectorView vector, bool transposeThis = false)
		{
			if (vector is Vector lhsDense)
			{
				return Multiply(lhsDense, transposeThis);
			}
			else
			{
				return base.Multiply(vector, transposeThis);
			}
		}

		/// <summary>
		/// Performs the matrix-vector multiplication: oper(this) * <paramref name="vector"/>.
		/// To multiply this * columnVector, set <paramref name="transposeThis"/> to false.
		/// To multiply rowVector * this, set <paramref name="transposeThis"/> to true.
		/// </summary>
		/// <param name="vector">A vector with <see cref="IIndexable1D.Length"/> being equal to <see cref="Order"/> of this 
		///     matrix.</param>
		/// <param name="transposeThis">If true, oper(this) = transpose(this). Otherwise oper(this) = this.</param>
		/// <exception cref="NonMatchingDimensionsException">Thrown if the <see cref="IIndexable1D.Length"/> of
		///     <paramref name="vector"/> is different than <see cref="Order"/> of this matrix.</exception>
		/// <returns>The result of the multiplication.</returns>
		public Vector Multiply(Vector vector, bool transposeThis = false)
		{
			//TODO: this performs redundant dimension checks
			var result = Vector.CreateZero(Order);
			MultiplyIntoResult(vector, result, transposeThis);
			return result;
		}

		public override void MultiplyIntoResult(IVectorView lhsVector, IVector rhsVector, bool transposeThis = false)
		{
			if ((lhsVector is Vector lhsDense) && (rhsVector is Vector rhsDense))
			{
				MultiplyIntoResult(lhsDense, rhsDense, transposeThis);
			}
			else
			{
				base.MultiplyIntoResult(lhsVector, rhsVector, transposeThis);
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
			TransposeMatrix transpose = transposeThis ? TransposeMatrix.NoTranspose : TransposeMatrix.Transpose; // row major
			Preconditions.CheckMultiplicationDimensions(Order, lhsVector.Length);
			Preconditions.CheckSystemSolutionDimensions(Order, rhsVector.Length);
			Array.Copy(lhsVector.RawData, rhsVector.RawData, Order);
			GlobalProvider.Blas.Dtpmv(
				StoredTriangle.Upper, transpose, DiagonalValues.NonUnit, Order, this.values, 0, rhsVector.RawData, 0, 1);
		}

		public override double Reduce(double identityValue, ProcessEntry processEntry, ProcessZeros processZeros, Finalize finalize)
			=> ReduceNonSymmetric(identityValue, processEntry, processZeros, finalize);

		/// <summary>
		/// Solves the linear equations system: this * result = <paramref name="rhs"/> by forward substitution. WARNING: this
		/// matrix must be invertible. No exception will be thrown if the matrix is singular.
		/// </summary>
		/// <param name="rhs">The right hand side vector of the linear system. Constraints: 
		///     <paramref name="rhs"/>.<see cref="Vector.Length"/> == this.<see cref="Order"/>.</param>
		public Vector SolveLinearSystem(Vector rhs, bool transposeThis = false)
		{
			Preconditions.CheckSystemSolutionDimensions(this, rhs);
			double[] result = rhs.CopyToArray();
			TransposeMatrix transposeBLAS = transposeThis ? TransposeMatrix.NoTranspose : TransposeMatrix.Transpose; // row major
			GlobalProvider.Blas.Dtpsv(
				StoredTriangle.Upper, transposeBLAS, DiagonalValues.NonUnit, Order, this.values, 0, result, 0, 1);
			return Vector.CreateFromArray(result, false);
		}

		public override IMatrix Transpose() => Transpose(true);

		/// <summary>
		/// Creates a new <see cref="TriangularUpper"/> matrix, that is transpose to this: result[i, j] = this[j, i]. The 
		/// internal array can be copied or shared with this <see cref="TriangularLower"/> matrix.
		/// </summary>
		/// <param name="copyInternalArray">If true, the internal array that stores the entries of this 
		///     <see cref="TriangularLower"/> instance will be copied and the new <see cref="TriangularUpper"/> instance 
		///     will have a reference to the copy, which is safer. If false, both the new matrix and this one will have  
		///     a reference to the same internal array, which is faster.</param>
		public TriangularUpper Transpose(bool copyInternalArray)
			=> TriangularUpper.CreateFromArray(Order, values, copyInternalArray); // trans(lower row major) = upper col major

		[MethodImpl(MethodImplOptions.AggressiveInlining)]
		private int FindIndex1D(int i, int j)
		{
			return j + ((i + 1) * i) / 2;
		}
	}
}
