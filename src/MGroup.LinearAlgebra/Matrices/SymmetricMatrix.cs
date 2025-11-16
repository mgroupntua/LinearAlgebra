//TODO: align data using mkl_malloc
namespace MGroup.LinearAlgebra.Matrices
{
	using System;
	using System.Runtime.CompilerServices;

	using MGroup.LinearAlgebra.Commons;
	using MGroup.LinearAlgebra.Eigensystems;
	using MGroup.LinearAlgebra.Exceptions;
	using MGroup.LinearAlgebra.Implementations;
	using MGroup.LinearAlgebra.Reduction;
	using MGroup.LinearAlgebra.Triangulation;
	using MGroup.LinearAlgebra.Vectors;

	using static MGroup.LinearAlgebra.LibrarySettings;

	/// <summary>
	/// Symmetric matrix. Only the upper triangle is stored in Packed format (only stores the n*(n+1)/2 non zeros) and column 
	/// major order. Uses LAPACK. Do not use this, since it is an experimantal class, which will probably be removed.
	/// </summary>
	[Serializable]
	public sealed class SymmetricMatrix : ValuesBackedMatrix<SymmetricMatrix>, ISymmetricMatrix,
		IEntrywiseOperableView2D<SymmetricMatrix, SymmetricMatrix>, IEntrywiseOperable2D<SymmetricMatrix>
	{
		/// <summary>
		/// Packed storage, column major order, upper triangle: 
		/// A[i,j] = data[i + j*(j+1)/2] for 0 &lt;= i &lt;= j &lt; n.
		/// </summary>
		private readonly double[] data;

		private SymmetricMatrix(double[] data, int order, DefiniteProperty definiteness)
		{
			this.data = data;
			this.Definiteness = definiteness;
			this.Order = order;
			this.NumRows = order;
			this.NumColumns = order;
		}

		public override MatrixSymmetry MatrixSymmetry => MatrixSymmetry.Symmetric;

		/// <summary>
		/// Used to query if the matrix is positive definite etc. Usually this is not known beforehand, which corresponds to
		/// <see cref="DefiniteProperty.Unknown"/>. Cholesky factorization reveals this and sets this property to
		/// <see cref="DefiniteProperty.PositiveDefinite"/> or <see cref="DefiniteProperty.Indefinite"/>. Mutating the matrix 
		/// will reset it to <see cref="DefiniteProperty.Unknown"/>. If the <see cref="SymmetricMatrix"/> is created from a 
		/// known matrix/array, the caller can assume responsibility for setting this property. WARNING: only set this propetry 
		/// if you are absolutely sure. 
		/// </summary>
		public DefiniteProperty Definiteness { get; set; }

		public override int NumRows { get; }

		public override int NumColumns { get; }

		/// <summary>
		/// The number of rows or columns of the matrix.
		/// </summary>
		public int Order { get; }

		/// <summary>
		/// The internal array that stores the entries of the upper triangle (packed storage format) in column major layout.
		/// It should only be used for passing the raw array to linear algebra libraries.
		/// </summary>
		public override double[] RawValues => data;

		public override double this[int rowIdx, int colIdx]
		{
			get
			{
				Preconditions.CheckIndices(this, rowIdx, colIdx);
				int index1D = (rowIdx <= colIdx) ? Find1DIndex(rowIdx, colIdx) : Find1DIndex(colIdx, rowIdx);
				return data[index1D];
			}

			set
			{
				Preconditions.CheckIndices(this, rowIdx, colIdx);
				int index1D = (rowIdx <= colIdx) ? Find1DIndex(rowIdx, colIdx) : Find1DIndex(colIdx, rowIdx);
				data[index1D] = value;
			}
		}

		/// <summary>
		/// Create a new <see cref="SymmetricMatrix"/> from the lower (subdiagonal) or upper (superdiagonal) portion of the 
		/// provided array. The array entries will be copied.
		/// </summary>
		/// <param name="array2D">A 2-dimensional containing the elements of the whole matrix. Its lengths in both dimensions 
		///     must be the same.</param>
		/// <param name="definiteness">If the caller knows that the matrix is positive definite, etc, he can set this property 
		///     during creation of the <see cref="SymmetricMatrix"/> object.</param>
		/// <returns></returns>
		public static SymmetricMatrix CreateFromArray(double[,] array2D,
			DefiniteProperty definiteness = DefiniteProperty.Unknown)
		{
			int numRows = array2D.GetLength(0);
			int numCols = array2D.GetLength(1);
			if (numRows != numCols)
			{
				string msg = string.Format("Provided array must have the same dimensions, but was ({0}x{1})", numRows, numCols);
				throw new NonMatchingDimensionsException(msg);
			}
			return new SymmetricMatrix(Conversions.Array2DToPackedUpperColMajor(array2D), numRows, definiteness);
		}

		/// <summary>
		/// Create a new <see cref="SymmetricMatrix"/> from a provided array. The array can be copied (for extra safety)
		/// or not (for extra performance).
		/// </summary>
		/// <param name="array1D">A 1-dimensional array containing the elements of the upper triangle of the matrix in column 
		///     major order.</param>
		/// <param name="order"> The order of the matrix. It must be positive and match the length of <see cref="array1D"/>. If a 
		///     value is provided, these will not be checked. If no value is provided, the order will be calculated from 
		///     <see cref="array1D"/> instead.</param>
		/// <param name="definiteness">If the caller knows that the matrix is positive definite etc., he can set this property 
		///     during creation of the <see cref="SymmetricMatrix"/> object.</param>
		/// <param name="copyArray">True to make a deep copy of <see cref="array1D"/>. 
		///     False (default) to use <see cref="array1D"/> as its internal storage.</param>
		public static SymmetricMatrix CreateFromPackedColumnMajorArray(double[] array1D, int order = 0,
			DefiniteProperty definiteness = DefiniteProperty.Unknown, bool copyArray = false)
		{
			int n = (order == 0) ? Conversions.PackedLengthToOrder(array1D.Length) : order;
			if (copyArray)
			{
				var clone = new double[array1D.Length];
				Array.Copy(array1D, clone, array1D.Length);
				return new SymmetricMatrix(clone, n, definiteness);
			}
			else return new SymmetricMatrix(array1D, n, definiteness);
		}

		/// <summary>
		/// Create a new <see cref="SymmetricMatrix"/> from a provided array. The array can be copied (for extra safety)
		/// or not (for extra performance).
		/// </summary>
		/// <param name="array1D">A 1-dimensional array containing the elements of the upper triangle of the matrix in row 
		///     major order.</param>
		/// <param name="order"> The order of the matrix. It must be positive and match the length of <see cref="array1D"/>. If a 
		///     value is provided, these will not be checked. If no value is provided, the order will be calculated from 
		///     <see cref="array1D"/> instead.</param>
		/// <param name="definiteness">If the caller knows that the matrix is positive definite etc., he can set this property 
		///     during creation of the <see cref="SymmetricMatrix"/> object.</param>
		public static SymmetricMatrix CreateFromPackedRowMajorArray(double[] array1D, int order = 0,
			DefiniteProperty definiteness = DefiniteProperty.Unknown)
		{
			int n = (order == 0) ? Conversions.PackedLengthToOrder(array1D.Length) : order;
			double[] columnMajor = Conversions.PackedUpperRowMajorPackedUpperColMajor(n, array1D);
			return new SymmetricMatrix(columnMajor, n, definiteness);
		}

		/// <summary>
		/// The caller is responsible for the original matrix being symmetric
		/// </summary>
		/// <param name="originalMatrix"></param>
		/// <returns></returns>
		public static SymmetricMatrix CreateFromMatrix(Matrix originalMatrix)
		{
			double[] data = Conversions.FullColMajorToPackedUpperColMajor(originalMatrix.NumColumns,
				originalMatrix.RawData);
			return new SymmetricMatrix(data, originalMatrix.NumColumns, DefiniteProperty.Unknown);
		}

		/// <summary>
		/// Create a new <see cref="SymmetricMatrix"/> with the specified order and all entries equal to 0.
		/// </summary> 
		/// <param name="order">The number of rows or columns of the matrix.</param>
		/// <returns></returns>
		public static SymmetricMatrix CreateZero(int order)
		{
			double[] data = new double[((order + 1) * order) / 2];
			//This matrix will be used as a canvas, thus we cannot infer that it is indefinite yet.
			return new SymmetricMatrix(data, order, DefiniteProperty.Unknown);
		}

		#region operators (use extension operators when they become available)
		public static SymmetricMatrix operator +(SymmetricMatrix matrix1, SymmetricMatrix matrix2)
			=> matrix1.AxpySameFormat(matrix2, 1.0);

		public static SymmetricMatrix operator -(SymmetricMatrix matrix1, SymmetricMatrix matrix2)
			=> matrix1.AxpySameFormat(matrix2, -1.0);

		public static SymmetricMatrix operator *(double scalar, SymmetricMatrix matrix)
			=> matrix.ScaleSameFormat(scalar);

		public static SymmetricMatrix operator *(SymmetricMatrix matrix, double scalar)
			=> matrix.ScaleSameFormat(scalar);

		public static IReadOnlyMatrix operator *(SymmetricMatrix matrixLeft, IReadOnlyMatrix matrixRight)
			=> matrixLeft.MultiplyRight(matrixRight, false, false);

		public static IReadOnlyMatrix operator *(IReadOnlyMatrix matrixLeft, SymmetricMatrix matrixRight)
			=> matrixRight.MultiplyLeft(matrixLeft, false, false);

		public static Vector operator *(SymmetricMatrix matrixLeft, Vector vectorRight)
			=> matrixLeft.Multiply(vectorRight);

		public static Vector operator *(Vector vectorLeft, SymmetricMatrix matrixRight)
			=> matrixRight.Multiply(vectorLeft);

		#endregion

		public override IMatrix Axpy(IReadOnlyMatrix otherMatrix, double otherCoefficient)
		{
			if (otherMatrix is SymmetricMatrix casted) return Axpy(casted, otherCoefficient);
			else return DoEntrywise(otherMatrix, (x1, x2) => x1 + otherCoefficient * x2); //TODO: optimize this
		}

		public override void AxpyIntoThis(SymmetricMatrix otherMatrix, double otherCoefficient)
		{
			base.AxpyIntoThis(otherMatrix, otherCoefficient);
			this.Definiteness = DefiniteProperty.Unknown;
		}

		/// <summary>
		/// Calculate the determinant of this matrix.
		/// If <see cref="Definiteness"/> != <see cref="DefiniteProperty.PositiveDefinite"/>, the calculation will be very 
		/// cumbersome and require an extra O(n^2) space: It will expand the packed storage format into a full data array, 
		/// factorize it using LU and then calculate the determinant. For positive definite matrices, it is more efficient,
		/// since it uses the Cholesky factorization directly.
		/// </summary>
		/// <returns></returns>
		public double CalcDeterminant()
		{
			if (Definiteness == DefiniteProperty.PositiveDefinite)
			{
				return FactorCholesky().CalcDeterminant();
			}
			else
			{
				return CopyToFullMatrix().FactorLU().CalcDeterminant(); //TODO: Find how to do it with Bunch-Kaufman
			}
		}

		/// <summary>
		/// Calculates the eigenvalues and eigenvectors of the symmetric matrix.
		/// </summary>
		/// <returns>The eigenvalues and eigenvectors of the matrix.</returns>
		public (Vector eigenvalues, Matrix eigenvectors) CalcEigensystem()
		{
			Preconditions.CheckSquare(this);
			double[] fullMatrix = Conversions.PackedUpperColMajorToFullSymmColMajor(data, Order);
			var eigensystem = SymmetricEigensystemFull.Create(Order, fullMatrix, true);
			return (eigensystem.EigenvaluesReal, eigensystem.EigenvectorsRight);
		}

		public override SymmetricMatrix CopyAsSameType(bool copyIndexingData = false)
		{
			double[] clone = new double[data.Length];
			Array.Copy(data, clone, data.Length);
			return new SymmetricMatrix(clone, Order, Definiteness);
		}

		/// <summary>
		/// Copy the entries of the matrix into a 2-dimensional array. The returned array has
		/// length(0) = <see cref="Order"/> and length(1) = <see cref="Order"/>.
		/// </summary>
		/// <returns>
		/// A new <see cref="double"/>[<see cref="Order"/>, <see cref="Order"/>] array with the entries of the matrix.
		/// </returns>
		public double[,] CopyToArray2D() => Conversions.PackedUpperColMajorToArray2DSymm(data, Order);

		public override Matrix CopyToFullMatrix()
		{
			double[] fullData = Conversions.PackedUpperColMajorToFullSymmColMajor(data, Order);
			return Matrix.CreateFromArray(fullData, Order, Order, false);
		}

		public override SymmetricMatrix CreateZeroMatrixSame()
		{
			var resultValues = new double[data.Length];
			return new SymmetricMatrix(resultValues, Order, DefiniteProperty.Unknown);
		}

		public SymmetricMatrix DoEntrywise(SymmetricMatrix matrix, Func<double, double, double> binaryOperation)
			=> DoEntrywiseSameFormat(matrix, binaryOperation);

		public override void DoEntrywiseIntoThis(IReadOnlyMatrix other, Func<double, double, double> binaryOperation)
		{
			if (other is SymmetricMatrix casted)
			{
				DoEntrywiseIntoThis(casted, binaryOperation);
			}
			else if (other is ISymmetricMatrix otherSYM)
			{
				Preconditions.CheckSameMatrixDimensions(this, other);
				for (int j = 0; j < NumColumns; ++j)
				{
					for (int i = 0; i <= j; ++i)
					{
						int index1D = Find1DIndex(i, j);
						this.data[index1D] = binaryOperation(this.data[index1D], other[i, j]);
					}
				}
			}
			else
			{
				throw new SymmetricPatternModifiedException("This operation is legal only if the other matrix is also symmetric.");
			}
		}

		public override void DoEntrywiseIntoThis(SymmetricMatrix otherMatrix, Func<double, double, double> binaryOperation)
		{
			base.DoEntrywiseIntoThis(otherMatrix, binaryOperation);
			this.Definiteness = DefiniteProperty.Unknown;
		}

		SymmetricMatrix IEntrywiseOperableView2D<SymmetricMatrix, SymmetricMatrix>.DoToAllEntries(
			Func<double, double> unaryOperation) => DoToAllEntriesSameFormat(unaryOperation);

		/// <summary>
		/// Calculates some factorization of the symmetric matrix.
		/// </summary>
		/// <param name="tryCholesky">True to apply the Cholesky factorization, before resorting to more expensive procedures.
		///     Since Cholesky has a complecity of O(1/3*n^3), this will further increase the computational cost for matrices 
		///     that are not positive definite. Therefore set it to false if you are sure that the matrix is not positive
		///     definite.</param>
		/// <returns></returns>
		public ITriangulation Factorize(bool tryCholesky = true)
		{
			// This should throw an exception if the posdef assumption is wrong.
			if (Definiteness == DefiniteProperty.PositiveDefinite) return FactorCholesky();

			if (tryCholesky)
			{
				try
				{
					return FactorCholesky();
				}
				catch (IndefiniteMatrixException) { }
			}

			return FactorBunchKaufman();
		}

		public BunchKaufmanFactorization FactorBunchKaufman()
		{
			throw new NotImplementedException();
		}

		/// <summary>
		/// Calculates the Cholesky factorization of the matrix. Will throw <see cref="IndefiniteMatrixException"/> if the 
		/// matrix is not positive definite.
		/// </summary>
		/// <returns></returns>
		public CholeskyPacked FactorCholesky()
		{
			var factor = CholeskyPacked.Factorize(Order, data);
			Definiteness = DefiniteProperty.PositiveDefinite; // An exception would have been thrown otherwise.
			return factor;
		}

		public override Vector GetColumn(int colIndex)
		{
			Preconditions.CheckIndexCol(this, colIndex);
			var columnVector = new double[Order];

			// Upper triangle and diagonal entries of the column are stored explicitly and contiguously
			int colOffset = (colIndex * (colIndex + 1)) / 2;
			Array.Copy(data, colOffset, columnVector, 0, colIndex + 1);

			// Lower triangle entries of the column can be found in the row with the same index
			for (int j = colIndex + 1; j < Order; ++j) columnVector[j] = data[colIndex + (j * (j + 1)) / 2];

			return Vector.CreateFromArray(columnVector);
		}

		public override Vector GetRow(int rowIndex) => GetColumn(rowIndex);

		public override bool HasSameFormat(SymmetricMatrix other) => this.Order == other.Order;

		public override void LinearCombinationIntoThis(double thisCoefficient, IReadOnlyMatrix otherMatrix, double otherCoefficient)
		{
			if (otherMatrix is SymmetricMatrix casted)
			{
				LinearCombinationIntoThis(thisCoefficient, casted, otherCoefficient);
			}
			else if (otherMatrix is ISymmetricMatrix otherSYM)
			{
				Preconditions.CheckSameMatrixDimensions(this, otherMatrix);
				for (int j = 0; j < NumColumns; ++j)
				{
					for (int i = 0; i <= j; ++i)
					{
						int index1D = Find1DIndex(i, j);
						this.data[index1D] = thisCoefficient * this.data[index1D] + otherCoefficient * otherMatrix[i, j];
					}
				}
			}
			else
			{
				throw new SymmetricPatternModifiedException(
					"This operation is legal only if the other matrix is also symmetric.");
			}
		}

		public override void LinearCombinationIntoThis(
			double thisCoefficient, SymmetricMatrix otherMatrix, double otherCoefficient)
		{
			base.LinearCombinationIntoThis(thisCoefficient, otherMatrix, otherCoefficient);
			this.Definiteness = DefiniteProperty.Unknown;
		}

		public override IVector Multiply(IReadOnlyVector vector, bool transposeThis = false)
		{
			if (vector is Vector lhsDense)
			{
				return Multiply(lhsDense);
			}
			else
			{
				return base.Multiply(vector, transposeThis);
			}
		}

		/// <summary>
		/// Matrix vector multiplication, with the vector on the right: matrix * vector.
		/// </summary>
		/// <param name="vector">A vector with length equal to <see cref="NumColumns"/>.</param>
		/// <returns></returns>
		public Vector Multiply(Vector vector)
		{
			//TODO: this performs redundant dimension checks
			var result = Vector.CreateZero(Order);
			MultiplyIntoResult(vector, result);
			return result;
		}

		public override void MultiplyIntoResult(IReadOnlyVector lhsVector, IVector rhsVector, bool transposeThis = false)
		{
			if ((lhsVector is Vector lhsDense) && (rhsVector is Vector rhsDense))
			{
				MultiplyIntoResult(lhsDense, rhsDense);
			}
			else
			{
				base.MultiplyIntoResult(lhsVector, rhsVector, transposeThis);
			}
		}

		/// <summary>
		/// Performs the matrix-vector multiplication: <paramref name="rhsVector"/> = this * <paramref name="vector"/>.
		/// The resulting vector will overwrite the entries of <paramref name="rhsVector"/>.
		/// </summary>
		/// <param name="lhsVector">
		/// The vector that will be multiplied by this matrix. It sits on the left hand side of the equation y = A * x.
		/// Constraints: <paramref name="lhsVector"/>.<see cref="IIndexable1D.Length"/> 
		/// == this.<see cref="IIndexable2D.NumColumns"/>.
		/// </param>
		/// <param name="rhsVector">
		/// The vector that will be overwritten by the result of the multiplication. It sits on the right hand side of the 
		/// equation y = A * x. Constraints: <paramref name="lhsVector"/>.<see cref="IIndexable1D.Length"/> 
		/// == this.<see cref="IIndexable2D.NumRows"/>.
		/// </param>
		/// <exception cref="NonMatchingDimensionsException">
		/// Thrown if the <see cref="IIndexable1D.Length"/> of <paramref name="lhsVector"/> or <paramref name="rhsVector"/> 
		/// violate the described contraints.
		/// </exception>
		public void MultiplyIntoResult(Vector lhsVector, Vector rhsVector)
		{
			Preconditions.CheckMultiplicationDimensions(this.NumColumns, lhsVector.Length);
			Preconditions.CheckSystemSolutionDimensions(this.NumRows, rhsVector.Length);
			GlobalProvider.Blas.Dspmv(StoredTriangle.Upper, Order,
				1.0, this.data, 0, lhsVector.RawData, 0, 1,
				0.0, rhsVector.RawData, 0, 1);
		}

		public override double Reduce(double identityValue, ProcessEntry processEntry, ProcessZeros processZeros, Finalize finalize)
		{
			double aggregator = identityValue;
			for (int j = 0; j < data.Length; ++j)
			{
				aggregator = processEntry(data[Find1DIndex(j, j)], aggregator);
				for (int i = 0; i < j; ++j)
				{
					// A[j,i] = A[i,j], but doubling it will not work for all reductions
					int idx1D = Find1DIndex(i, j);
					aggregator = processEntry(data[idx1D], aggregator);
					aggregator = processEntry(data[idx1D], aggregator);
				}
			}
			// no zeros implied
			return finalize(aggregator);
		}

		public override IMatrix Transpose() => Transpose(true);

		public SymmetricMatrix Transpose(bool copyInternalArray) 
			=> CreateFromPackedColumnMajorArray(data, Order, Definiteness, copyInternalArray);

		[MethodImpl(MethodImplOptions.AggressiveInlining)]
		internal int Find1DIndex(int i, int j) => i + (j * (j + 1)) / 2;
	}
}
