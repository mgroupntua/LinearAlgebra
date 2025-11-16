//TODO: Perhaps I should use row major for lower triangular, upper triangular or both.
//TODO: Perhaps I should have an abstract class that handles everything except the lower/upper specific stuff and concrete
//  private classes Lower, Upper. The indexer would be faster.
//TODO: align data using mkl_malloc
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
	/// Upper triangular square matrix in column major Packed storage format (only stores the n*(n+1)/2 non zeros). Uses LAPACK.  
	/// For the more information about the layout, see 
	/// <see cref="https://software.intel.com/en-us/mkl-developer-reference-c-matrix-storage-schemes-for-lapack-routines."/>
	/// </summary>
	public sealed class TriangularUpper : ValuesBackedMatrix<TriangularUpper>
	{
		/// <summary>
		/// Packed storage, column major order: U[i, j] = data[i + j*(2*n-j-1)/2] for 0 &lt;= j &lt;= i &lt; n.
		/// Although not used here, lower triangular would be: L[i, j] = data[i + j*(j+1)/2] for 0 &lt;= i &lt;= j &lt; n.
		/// </summary>
		private readonly double[] values;

		private TriangularUpper(double[] data, int order)
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
		/// The internal array that stores the entries of the upper triangle (packed storage format) in column major layout. 
		/// It should only be used for passing the raw array to linear algebra libraries.
		/// </summary>
		public override double[] RawValues => values;

		public override double this[int rowIdx, int colIdx]
		{
			get
			{
				Preconditions.CheckIndices(this, rowIdx, colIdx);
				return (colIdx >= rowIdx) ? values[FindIndex1D(rowIdx, colIdx)] : 0.0;
			}

			set
			{
				Preconditions.CheckIndices(this, rowIdx, colIdx);
				if (colIdx >= rowIdx)
				{
					values[FindIndex1D(rowIdx, colIdx)] = value;
				}
				else
				{
					throw new SparsityPatternModifiedException(
						$"Cannot change the subdiagonal entry A[{rowIdx}, {colIdx}] = 0.");
				}
			}
		}

		/// <summary>
		/// Initializes a new instance of <see cref="TriangularUpper"/> by copying the upper triangle of the provided 2D array.
		/// </summary>
		/// <param name="array2D">A 2-dimensional array containing the elements of the matrix. Constraints: 
		///     <paramref name="array2D"/>.GetLength(0) == <paramref name="array2D"/>.GetLength(1).</param>
		/// <exception cref="NonMatchingDimensionsException">Thrown if <paramref name="array2D"/>.GetLength(0) != 
		///     <paramref name="array2D"/>.GetLength(1).</exception>
		public static TriangularUpper CreateFromArray(double[,] array2D)
		{
			int numRows = array2D.GetLength(0);
			int numCols = array2D.GetLength(1);
			if (numRows != numCols)
			{
				string msg = $"Provided array must have the same dimensions, but was ({numRows}x{numCols})";
				throw new NonMatchingDimensionsException(msg);
			}
			return new TriangularUpper(Conversions.Array2DToPackedUpperColMajor(array2D), numRows);
		}

		/// <summary>
		/// Initializes a new instance of <see cref="TriangularUpper"/> with <paramref name="array1D"/> or a clone as its 
		/// internal array.
		/// </summary>
		/// <param name="order">The number of rows/columns of the new square matrix.</param>
		/// <param name="array1D">An 1-dimensional array containing the elements of the upper triangle of the matrix in column 
		///     major order.</param>
		/// <param name="copyArray">If true, <paramref name="array1D"/> will be copied and the new <see cref="TriangularUpper"/>  
		///     instance will have a reference to the copy, which is safer. If false, the new matrix will have a reference to 
		///     <paramref name="array1D"/> itself, which is faster.</param>
		public static TriangularUpper CreateFromArray(int order, double[] array1D, bool copyArray = false)
		{
			if (copyArray)
			{
				var clone = new double[array1D.Length];
				Array.Copy(array1D, clone, array1D.Length);
				return new TriangularUpper(clone, order);
			}
			else return new TriangularUpper(array1D, order);
		}

		/// <summary>
		/// Initializes a new instance of <see cref="TriangularUpper"/> with all entries being equal to 0. Only the upper 
		/// triangle zero entries are explictily stored though.
		/// </summary>
		/// <param name="order">The number of rows/columns of the new square matrix.</param>
		public static TriangularUpper CreateZero(int order)
		{
			var data = new double[(order * (order + 1)) / 2];
			return new TriangularUpper(data, order);
		}

		/// <summary>
		/// Calculates the determinant of this matrix.
		/// </summary>
		public double CalcDeterminant()
		{
			// TODO: Find more effienct formulas for the diagonal accesses.
			double det = 1.0;
			for (int i = 0; i < Order; ++i) det *= values[FindIndex1D(i, i)];
			return det;
		}

		public override TriangularUpper CopyAsSameType(bool copyIndexingData = false)
		{
			var clone = new double[this.values.Length];
			Array.Copy(this.values, clone, this.values.Length);
			return new TriangularUpper(clone, Order);
		}

		/// <summary>
		/// Copies the entries of the matrix into a 2-dimensional array. The returned array has
		/// length(0) = length(1) = <see cref="Order"/>.
		/// </summary>
		public double[,] CopyToArray2D() => Conversions.PackedUpperColMajorToArray2D(values);

		public override Matrix CopyToFullMatrix()
		{
			// TODO: This won't work if the implementation of Matrix changes
			double[] fullArray = Conversions.PackedUpperColMajorToFullColMajor(values, Order);
			return Matrix.CreateFromArray(fullArray, Order, Order, false);
		}

		public override TriangularUpper CreateZeroMatrixSame()
		{
			var resultValues = new double[values.Length];
			return new TriangularUpper(resultValues, Order);
		}

		public override Vector GetColumn(int colIndex)
		{
			Preconditions.CheckIndexCol(this, colIndex);
			var columnVector = new double[NumRows];
			int numNonZerosCol = colIndex + 1;
			int colOffset = (colIndex * (colIndex + 1)) / 2;
			Array.Copy(values, colOffset, columnVector, 0, numNonZerosCol);
			return Vector.CreateFromArray(columnVector);
		}

		public override Vector GetRow(int rowIndex)
		{
			Preconditions.CheckIndexRow(this, rowIndex);
			var rowVector = new double[NumColumns];
			for (int j = rowIndex; j < NumRows; ++j) rowVector[j] = values[rowIndex + (j * (j + 1)) / 2];
			return Vector.CreateFromArray(rowVector);
		}

		public override bool HasSameFormat(TriangularUpper other) => this.Order == other.Order;

		public override IVector Multiply(IReadOnlyVector vector, bool transposeThis = false)
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
		public Vector Multiply(Vector vector, bool transposeThis = false)
		{
			//TODO: this performs redundant dimension checks
			var result = Vector.CreateZero(Order);
			MultiplyIntoResult(vector, result, transposeThis);
			return result;
		}

		public override void MultiplyIntoResult(IReadOnlyVector lhsVector, IVector rhsVector, bool transposeThis = false)
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
			TransposeMatrix transpose = transposeThis ? TransposeMatrix.Transpose : TransposeMatrix.NoTranspose;
			Preconditions.CheckMultiplicationDimensions(Order, lhsVector.Length);
			Preconditions.CheckSystemSolutionDimensions(Order, rhsVector.Length);
			Array.Copy(lhsVector.RawData, rhsVector.RawData, Order);
			GlobalProvider.Blas.Dtpmv(StoredTriangle.Upper, transpose, DiagonalValues.NonUnit, Order,
				this.values, 0, rhsVector.RawData, 0, 1);
		}

		public override double Reduce(double identityValue, ProcessEntry processEntry, ProcessZeros processZeros, Finalize finalize)
			=> ReduceNonSymmetric(identityValue, processEntry, processZeros, finalize);

		/// <summary>
		/// Solves the linear equations system: this * result = <paramref name="rhs"/> by back substitution. WARNING: this
		/// matrix must be invertible. No exception will be thrown if the matrix is singular.
		/// </summary>
		/// <param name="rhs">The right hand side vector of the linear system. Constraints: 
		///     <paramref name="rhs"/>.<see cref="Vector.Length"/> == this.<see cref="Order"/>.</param>
		public Vector SolveLinearSystem(Vector rhs, bool transposeThis = false)
		{
			Preconditions.CheckSystemSolutionDimensions(this, rhs);
			double[] result = rhs.CopyToArray();
			TransposeMatrix transposeBlas = transposeThis ? TransposeMatrix.Transpose : TransposeMatrix.NoTranspose;
			GlobalProvider.Blas.Dtpsv(StoredTriangle.Upper, transposeBlas, DiagonalValues.NonUnit, Order,
				this.values, 0, result, 0, 1);
			return Vector.CreateFromArray(result, false);
		}

		public override IMatrix Transpose() => Transpose(true);

		/// <summary>
		/// Creates a new <see cref="TriangularLower"/> matrix, that is transpose to this: result[i, j] = this[j, i]. The  
		/// internal array can be copied or shared with this <see cref="TriangularUpper"/> matrix.
		/// </summary>
		/// <param name="copyInternalArray">If true, the internal array that stores the entries of this 
		///     <see cref="TriangularUpper"/> instance will be copied and the new <see cref="TriangularLower"/> instance 
		///     will have a reference to the copy, which is safer. If false, both the new matrix and this one will have  
		///     a reference to the same internal array, which is faster.</param>
		public TriangularLower Transpose(bool copyInternalArray)
			=> TriangularLower.CreateFromArray(Order, values, copyInternalArray); // trans(upper col major) = lower row major


		[MethodImpl(MethodImplOptions.AggressiveInlining)]
		private int FindIndex1D(int i, int j)
		{
			return i + ((j + 1) * j) / 2;
		}
	}
}
