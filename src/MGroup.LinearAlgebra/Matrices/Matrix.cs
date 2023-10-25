using System;
using System.Collections.Generic;
using MGroup.LinearAlgebra.Commons;
using MGroup.LinearAlgebra.Exceptions;
using MGroup.LinearAlgebra.Triangulation;
using MGroup.LinearAlgebra.Providers;
using MGroup.LinearAlgebra.Reduction;
using MGroup.LinearAlgebra.Vectors;
using static MGroup.LinearAlgebra.LibrarySettings;
using MGroup.LinearAlgebra.Orthogonalization;
using MGroup.LinearAlgebra.Eigensystems;

//TODO: align data using mkl_malloc
//TODO: add inplace option for factorizations and leave all subsequent operations (determinant, system solution, etc.) to them
//TODO: remove legacy matrix conversions
//TODO: SetSubrow, SetSubcolumn, SetSubmatrix only need to check the stricter upper bounds.
//TODO: Se https://software.intel.com/en-us/mkl-developer-reference-c-lapmr, 
//      https://software.intel.com/en-us/mkl-developer-reference-c-laswp, 
//      https://software.intel.com/en-us/mkl-developer-reference-c-syswapr for reordering
namespace MGroup.LinearAlgebra.Matrices
{
	/// <summary>
	/// General purpose matrix class. All entries are stored in an 1D column major array. Uses LAPACK for most operations. 
	/// Authors: Serafeim Bakalakos
	/// </summary>
	[Serializable]
	public class Matrix : IMatrix, ISliceable2D, IEntrywiseOperableView2D<Matrix, Matrix>, IEntrywiseOperable2D<Matrix>
	{
		private double[] data;
		private bool isOverwritten = false;

		private Matrix(double[] data, int numRows, int numColumns)
		{
			this.data = data;
			this.NumRows = numRows;
			this.NumColumns = numColumns;
		}

		/// <summary>
		/// Returns true if <see cref="NumRows"/> == <see cref="NumColumns"/>.
		/// </summary>
		public bool IsSquare { get { return NumRows == NumColumns; } }

		/// <summary>
		/// See <see cref="IIndexable2D.MatrixSymmetry"/>.
		/// </summary>
		public MatrixSymmetry MatrixSymmetry { get; set; }

		/// <summary>
		/// See <see cref="IIndexable2D.MatrixSymmetry"/>.
		/// </summary>
		MatrixSymmetry IIndexable2D.MatrixSymmetry => this.MatrixSymmetry;

		/// <summary>
		/// The number of columns of the matrix. 
		/// </summary>
		public int NumColumns { get; }

		/// <summary>
		/// The number of non-zero and explicitly stored zero entries, which is the number of all entries in this 
		/// <see cref="Matrix"/>.
		/// </summary>
		public int NumNonZeros { get { return NumRows * NumColumns; } }

		/// <summary>
		/// The number of rows of the matrix.
		/// </summary>
		public int NumRows { get; }

		/// <summary>
		/// The internal array that stores the entries of the matrix in column major layout.
		/// It should only be used for passing the raw array to linear algebra libraries.
		/// </summary>
		public double[] RawData { get { return data; } }

		/// <summary>
		/// See <see cref="IIndexable2D.this[int, int]"/>.
		/// </summary>
		/// <remarks>
		/// Also note that it may be possible to pass in <paramref name="rowIdx"/> &gt;= <see cref="IIndexable2D.NumRows"/> or
		/// <paramref name="rowIdx"/> &lt; 0, without throwing <see cref="IndexOutOfRangeException"/>, since the indices are not  
		/// checked explicitly. The constraints on <paramref name="colIdx"/> described in the interfaces will correctly throw
		/// <see cref="IndexOutOfRangeException"/> if violated.
		/// </remarks>
		public double this[int rowIdx, int colIdx] //TODO: Should I add bound checking?
		{
			get { return data[colIdx * NumRows + rowIdx]; }
			set { data[colIdx * NumRows + rowIdx] = value; }
		}

		/// <summary>
		/// Initializes a new instance of <see cref="Matrix"/> by copying the entries of <paramref name="array2D"/>.
		/// </summary>
		/// <param name="array2D">A 2-dimensional array containing the entries of the matrix. It will be copied.</param>
		public static Matrix CreateFromArray(double[,] array2D)
		{
			int numRows = array2D.GetLength(0);
			int numCols = array2D.GetLength(1); 
			return new Matrix(Conversions.Array2DToFullColMajor(array2D), numRows, numCols);
		}

		/// <summary>
		/// Initializes a new instance of <see cref="Matrix"/> with <paramref name="array1D"/> or a clone as its internal array.
		/// </summary>
		/// <param name="array1D">A 1-dimensional array containing the elements of the matrix in column major order. Its length 
		///     must be equal to <see cref="numRows"/> * <see cref="NumColumns"/>. It will not be checked.</param>
		/// <param name="numRows">The number of rows of the new matrix.</param>
		/// <param name="numColumns">The number of columns of the new matrix.</param>
		/// <param name="copyArray">If true, <paramref name="array1D"/> will be copied and the new <see cref="Matrix"/> instance 
		///     will have a reference to the copy, which is safer. If false, the new matrix will have a reference to 
		///     <paramref name="array1D"/> itself, which is faster.</param>
		public static Matrix CreateFromArray(double[] array1D, int numRows, int numColumns, bool copyArray = false)
		{
			if (copyArray)
			{
				var clone = new double[array1D.Length];
				Array.Copy(array1D, clone, clone.Length);
				return new Matrix(clone, numRows, numColumns);
			}
			else
			{
				return new Matrix(array1D, numRows, numColumns);
			}
		}

		/// <summary>
		/// Initializes a new instance of <see cref="Matrix"/> with the entries of <paramref name="diagonal"/> on its main 
		/// diagonal.
		/// </summary>
		/// <param name="numRows">The number of rows of the new matrix.</param>
		/// <param name="numColumns">The number of columns of the new matrix.</param>
		/// <param name="diagonal">
		/// An array containing the entries of the main diagonal of the new matrix. Its length must be equal to 
		/// min(<paramref name="numRows"/>, <paramref name="numColumns"/>).
		/// </param>
		/// <returns></returns>
		public static Matrix CreateFromDiagonal(int numRows, int numColumns, double[] diagonal)
		{
			Preconditions.CheckVectorDimensions(Math.Min(numRows, numColumns), diagonal.Length);
			double[] data = new double[numRows * numColumns];
			ArrayColMajor.DiagonalSet(numRows, numColumns, data, diagonal);
			return new Matrix(data, numRows, numColumns);
		}

		/// <summary>
		/// Initializes a new instance of <see cref="Matrix"/> that is equal to the identity matrix, namely a square matrix with 
		/// non-diagonal entries being equal to 0 and diagonal entries being equal to 1.
		/// </summary>
		/// <param name="order">The number of rows/columns of the identity matrix.</param>
		public static Matrix CreateIdentity(int order)
		{
			double[] data = new double[order * order];
			for (int j = 0; j < order; ++j) data[j * order + j] = 1.0;
			return new Matrix(data, order, order);
		}

		/// <summary>
		/// Initializes a new instance of <see cref="Matrix"/> with all entries being equal to <paramref name="value"/>.
		/// </summary>
		/// <param name="numRows">The number of rows of the new matrix.</param>
		/// <param name="numColumns">The number of columns of the new matrix.</param>
		/// <param name="value">The value that all entries of the new matrix will be initialized to.</param>
		public static Matrix CreateWithValue(int numRows, int numColumns, double value)
		{
			double[] data = new double[numRows * numColumns];
			for (int i = 0; i < data.Length; ++i) data[i] = value;
			return new Matrix(data, numRows, numColumns);
		}

		/// <summary>
		/// Initializes a new instance of <see cref="Matrix"/> with all entries being equal to 0.
		/// </summary> 
		/// <param name="numRows">The number of rows of the new matrix.</param>
		/// <param name="numColumns">The number of rows of the new matrix.</param>
		/// <returns></returns>
		public static Matrix CreateZero(int numRows, int numColumns)
		{
			double[] data = new double[numRows * numColumns];
			return new Matrix(data, numRows, numColumns);
		}

		#region operators (use extension operators when they become available)
		/// <summary>
		/// Performs the operation: result[i, j] = <paramref name="matrix1"/>[i, j] + <paramref name="matrix2"/>[i, j], 
		/// for 0 &lt;= i &lt; <see cref="NumRows"/>, 0 &lt;= j &lt; <see cref="NumColumns"/>.
		/// The resulting entries are written to a new <see cref="Matrix"/> instance.
		/// </summary>
		/// <param name="matrix1">The first <see cref="Matrix"/> operand. It must have as many rows and columns as 
		///     <paramref name="matrix2"/>.</param>
		/// <param name="matrix2">The second <see cref="Matrix"/> operand. It must have as many rows and columns as 
		///     <paramref name="matrix1"/>.</param>
		/// <exception cref="NonMatchingDimensionsException">Thrown if <paramref name="matrix1"/> and <paramref name="matrix2"/>
		///     have a different number of <see cref="NumRows"/> or <see cref="NumColumns"/>.</exception>
		public static Matrix operator +(Matrix matrix1, Matrix matrix2) => matrix1.Axpy(matrix2, 1.0);

		/// <summary>
		/// Performs the operation: result[i, j] = <paramref name="matrix1"/>[i, j] - <paramref name="matrix2"/>[i, j], 
		/// for 0 &lt;= i &lt; <see cref="NumRows"/>, 0 &lt;= j &lt; <see cref="NumColumns"/>.
		/// The resulting entries are written to a new <see cref="Matrix"/> instance.
		/// </summary>
		/// <param name="matrix1">The first <see cref="Matrix"/> operand. It must have as many rows and columns as 
		///     <paramref name="matrix2"/>.</param>
		/// <param name="matrix2">The second <see cref="Matrix"/> operand. It must have as many rows and columns as 
		///     <paramref name="matrix1"/>.</param>
		/// <exception cref="NonMatchingDimensionsException">Thrown if <paramref name="matrix1"/> and <paramref name="matrix2"/>
		///     have a different number of <see cref="NumRows"/> or <see cref="NumColumns"/>.</exception>
		public static Matrix operator -(Matrix matrix1, Matrix matrix2) => matrix1.Axpy(matrix2, -1.0);

		/// <summary>
		/// Performs the operation: result[i, j] = <paramref name="scalar"/> * <paramref name="matrix1"/>[i, j],
		/// for 0 &lt;= i &lt; <see cref="NumRows"/>, 0 &lt;= j &lt; <see cref="NumColumns"/>.
		/// The resulting entries are written to a new <see cref="Matrix"/> instance.
		/// </summary>
		/// <param name="scalar">The scalar value that will be multiplied with all vector entries.</param>
		/// <param name="matrix">The matrix to multiply.</param>
		public static Matrix operator *(double scalar, Matrix matrix) => matrix.Scale(scalar);

		/// <summary>
		/// Performs the operation: result[i, j] = <paramref name="scalar"/> * <paramref name="matrix1"/>[i, j],
		/// for 0 &lt;= i &lt; <see cref="NumRows"/>, 0 &lt;= j &lt; <see cref="NumColumns"/>.
		/// The resulting entries are written to a new <see cref="Matrix"/> instance.
		/// </summary>
		/// <param name="matrix">The matrix to multiply.</param>
		/// <param name="scalar">The scalar value that will be multiplied with all vector entries.</param>
		public static Matrix operator *(Matrix matrix, double scalar)=> matrix.Scale(scalar);

		/// <summary>
		/// Performs the matrix-matrix multiplication: result = <paramref name="matrixLeft"/> * <paramref name="matrixRight"/>.
		/// If <paramref name="matrixLeft"/> is m1-by-n1 and <paramref name="matrixRight"/> is m2-by-n2, then n1 must be equal to
		/// m2. The result will be an m1-by-n2 matrix, written to a new <see cref="Matrix"/> instance.
		/// </summary>
		/// <param name="matrixLeft">The <see cref="Matrix"/> operand on the left.</param>
		/// <param name="matrixRight">The <see cref="Matrix"/> operand on the right.</param>
		/// <exception cref="NonMatchingDimensionsException">Thrown if <paramref name="matrixLeft"/>.<see cref="NumColumns"/> is 
		///     different than <paramref name="matrixRight"/>.<see cref="NumRows"/>.</exception>
		public static Matrix operator *(Matrix matrixLeft, Matrix matrixRight)
			=> matrixLeft.MultiplyRight(matrixRight, false, false);

		/// <summary>
		/// Performs the matrix-vector multiplication: result = <paramref name="matrixLeft"/> * <paramref name="vectorRight"/>.
		/// If <paramref name="matrixLeft"/> is m1-by-n1 and <paramref name="vectorRight"/> has length = n2, then n1 must be 
		/// equal to n2. The result will be a vector with length = m1, written to a new <see cref="Vector"/> instance.
		/// </summary>
		/// <param name="matrixLeft">The <see cref="Matrix"/> operand on the left.</param>
		/// <param name="vectorRight">The <see cref="Vector"/> operand on the right. It can be considered as a column 
		///     vector.</param>
		/// <exception cref="NonMatchingDimensionsException">Thrown if <paramref name="matrixLeft"/>.<see cref="NumColumns"/> is 
		///     different than <paramref name="vectorRight"/>.<see cref="Vector.Length"/>.</exception>
		public static Vector operator *(Matrix matrixLeft, Vector vectorRight)
			=> matrixLeft.Multiply(vectorRight, false);

		/// <summary>
		/// Performs the matrix-vector multiplication: result = <paramref name="vectorLeft"/> * <paramref name="matrixRight"/>.
		/// If <paramref name="matrixRight"/> is m1-by-n1 and <paramref name="vectorLeft"/> has length = n2, then m1 must be 
		/// equal to n2. The result will be a vector with length = n1, written to a new <see cref="Vector"/> instance.
		/// </summary>
		/// <param name="vectorLeft">The <see cref="Vector"/> operand on the left. It can be considered as a row vector.</param>
		/// <param name="matrixRight">The <see cref="Matrix"/> operand on the right.</param>
		/// <exception cref="NonMatchingDimensionsException">Thrown if <paramref name="matrixRight"/>.<see cref="NumRows"/> is 
		///     different than <paramref name="vectorLeft"/>.<see cref="Vector.Length"/>.</exception>
		public static Vector operator *(Vector vectorLeft, Matrix matrixRight)
			=> matrixRight.Multiply(vectorLeft, true);
		#endregion

		/// <summary>
		/// Creates a new <see cref="Matrix"/> that contains all rows of this <see cref="Matrix"/> instance, followed by all rows 
		/// of <paramref name="matrix"/>. If this is m1-by-n1 and <paramref name="matrix"/> is m2-by-n2, then n2 must be equal 
		/// to n1 and the resulting matrix will be (m1+m2)-by-n1.
		/// </summary>
		/// <param name="matrix">The matrix whose rows will be appended after all rows of this <see cref="Matrix"/> 
		///     instance.</param>
		public Matrix AppendBottom(Matrix matrix)
		{
			Preconditions.CheckSameColDimension(this, matrix);
			double[] result = ArrayColMajor.JoinVertically(this.NumRows, this.NumColumns, this.data,
				matrix.NumRows, matrix.NumColumns, matrix.data);
			return new Matrix(result, this.NumRows + matrix.NumRows, NumColumns);
		}

		/// <summary>
		/// Creates a new <see cref="Matrix"/> that contains all columns of this <see cref="Matrix"/> instance, followed by all 
		/// columns of <paramref name="matrix"/>. If this is m1-by-n1 and <paramref name="matrix"/> is m2-by-n2, then m2 must be 
		/// equal to mn1 and the resulting matrix will be m1-by-(n1+n2).
		/// </summary>
		/// <param name="matrix">The matrix whose columns will be appended after all columns of this <see cref="Matrix"/> 
		///     instance.</param>
		public Matrix AppendRight(Matrix matrix)
		{
			Preconditions.CheckSameRowDimension(this, matrix);
			double[] result = ArrayColMajor.JoinHorizontally(this.NumRows, this.NumColumns, this.data,
				matrix.NumRows, matrix.NumColumns, matrix.data);
			return new Matrix(result, NumRows, this.NumColumns + matrix.NumColumns);
		}

		/// <summary>
		/// See <see cref="IMatrixView.Axpy(IMatrixView, double)"/>.
		/// </summary>
		public IMatrix Axpy(IMatrixView otherMatrix, double otherCoefficient)
		{
			if (otherMatrix is Matrix dense) return Axpy(dense, otherCoefficient);
			else return otherMatrix.LinearCombination(otherCoefficient, this, 1.0); // To avoid accessing zero entries
		}

		/// <summary>
		/// Performs the following operation for 0 &lt;= i &lt; <see cref="NumRows"/>, 0 &lt;= j &lt; <see cref="NumColumns"/>:
		/// result[i, j] = <paramref name="otherCoefficient"/> * <paramref name="otherMatrix"/>[i, j] + this[i, j]. 
		/// The resulting matrix is written to a new <see cref="Matrix"/> and then returned.
		/// </summary>
		/// <param name="otherMatrix">A matrix with the same <see cref="NumRows"/> and <see cref="NumColumns"/> as this 
		///     <see cref="Matrix"/> instance.</param>
		/// <param name="otherCoefficient">A scalar that multiplies each entry of <paramref name="otherMatrix"/>.</param>
		/// <exception cref="NonMatchingDimensionsException">Thrown if <paramref name="otherMatrix"/> has different 
		///     <see cref="NumRows"/> or <see cref="NumColumns"/> than this instance.</exception>
		public Matrix Axpy(Matrix otherMatrix, double otherCoefficient)
		{
			Preconditions.CheckSameMatrixDimensions(this, otherMatrix);
			//TODO: Perhaps this should be done using mkl_malloc and BLAS copy. 
			double[] result = new double[data.Length];
			Array.Copy(this.data, result, data.Length);
			Blas.Daxpy(data.Length, otherCoefficient, otherMatrix.data, 0, 1, result, 0, 1);
			return new Matrix(result, NumRows, NumColumns);
		}

		/// <summary>
		/// See <see cref="IMatrix.AxpyIntoThis(IMatrixView, double)"/>.
		/// </summary>
		public void AxpyIntoThis(IMatrixView otherMatrix, double otherCoefficient)
		{
			if (otherMatrix is Matrix dense) AxpyIntoThis(dense, otherCoefficient);
			else
			{
				Preconditions.CheckSameMatrixDimensions(this, otherMatrix);
				for (int j = 0; j < NumColumns; ++j)
				{
					for (int i = 0; i < NumRows; ++i)
					{
						this.data[j * NumRows + i] += otherCoefficient * otherMatrix[i, j];
					}
				}
			}
		}

		/// <summary>
		/// Performs the following operation for 0 &lt;= i &lt; <see cref="NumRows"/>, 0 &lt;= j &lt; <see cref="NumColumns"/>:
		/// this[i, j] = <paramref name="otherCoefficient"/> * <paramref name="otherMatrix"/>[i, j] + this[i, j]. 
		/// The resulting matrix overwrites the entries of this <see cref="Matrix"/> instance.
		/// </summary>
		/// <param name="otherMatrix">A matrix with the same <see cref="NumRows"/> and <see cref="NumColumns"/> as this 
		///     <see cref="Matrix"/> instance.</param>
		/// <param name="otherCoefficient">A scalar that multiplies each entry of <paramref name="otherMatrix"/>.</param>
		/// <exception cref="NonMatchingDimensionsException">Thrown if <paramref name="otherMatrix"/> has different 
		///     <see cref="NumRows"/> or <see cref="NumColumns"/> than this instance.</exception>
		public void AxpyIntoThis(Matrix otherMatrix, double otherCoefficient)
		{
			Preconditions.CheckSameMatrixDimensions(this, otherMatrix);
			Blas.Daxpy(data.Length, otherCoefficient, otherMatrix.data, 0, 1, this.data, 0, 1);
		}

		/// <summary>
		/// Calculates the determinant of this matrix, which must be square. If the inverse matrix is also needed, use
		/// <see cref="InvertAndDeterminant"/> instead.
		/// </summary>
		/// <exception cref="NonMatchingDimensionsException">Thrown if this matrix is not square.</exception>
		public double CalcDeterminant()
		{
			if ((NumRows == 2) && (NumColumns == 2))
			{
				return AnalyticFormulas.Matrix2x2ColMajorDeterminant(data);
			}
			else if ((NumRows == 3) && (NumColumns == 3))
			{
				return AnalyticFormulas.Matrix3x3ColMajorDeterminant(data);
			}
			else return FactorLU().CalcDeterminant();
		}

		/// <summary>
		/// Calculates the eigenvalues and eigenvectors of the matrix. The matrix can be symmetric or not, but it must be square.
		/// </summary>
		/// <returns>The eigenvalues and eigenvectors of the matrix.</returns>
		public (Vector eigenvaluesReal, Vector eigenvaluesImaginary, Matrix eigenvectorsRight, Matrix eigenvectorsLeft)
			CalcEigensystemNonSymmetric()
		{
			Preconditions.CheckSquare(this);
			double[] clone = CopyInternalData();
			var eigensystem = NonSymmetricEigensystemFull.Create(NumRows, clone, true, true);
			return (eigensystem.EigenvaluesReal, eigensystem.EigenvaluesImaginary,
				eigensystem.EigenvectorsRight, eigensystem.EigenvectorsLeft);
		}

		/// <summary>
		/// Calculates the eigenvalues and eigenvectors of the matrix. The matrix must be symmetric for this to work correctly.
		/// </summary>
		/// <returns>The eigenvalues and eigenvectors of the matrix.</returns>
		public (Vector eigenvalues, Matrix eigenvectors) CalcEigensystemSymmetric()
		{
			Preconditions.CheckSquare(this);
			double[] clone = CopyInternalData();
			var eigensystem = SymmetricEigensystemFull.Create(NumRows, clone, true);
			return (eigensystem.EigenvaluesReal, eigensystem.EigenvectorsRight);
		}

		/// <summary>
		/// See <see cref="IMatrix.Clear"/>.
		/// </summary>
		public void Clear() => Array.Clear(data, 0, data.Length);

		/// <summary>
		/// See <see cref="IMatrixView.Copy(bool)"/>.
		/// </summary>
		IMatrix IMatrixView.Copy(bool copyIndexingData) => Copy();

		/// <summary>
		/// Initializes a new instance of <see cref="Matrix"/> by copying the entries of this instance.
		/// </summary>
		public Matrix Copy()
		{
			//TODO: Perhaps this should use BLAS. 
			double[] clone = new double[data.Length];
			Array.Copy(data, clone, data.Length);
			return new Matrix(clone, NumRows, NumColumns);
		}

		/// <summary>
		/// Copies the entries of the matrix into a 2-dimensional array. The returned array has length(0) = <see cref="NumRows"/> 
		/// and length(1) = <see cref="NumColumns"/>. 
		/// </summary>
		public double[,] CopyToArray2D()
		{
			return Conversions.FullColMajorToArray2D(data, NumRows, NumColumns);
		}

		/// See <see cref="IMatrixView.CopyToFullMatrix()"/>
		/// </summary>
		public Matrix CopyToFullMatrix() => Copy();

		/// <summary>
		/// <summary>
		/// See <see cref="IMatrixView.DoEntrywise(IMatrixView, Func{double, double, double})"/>.
		/// </summary>
		public IMatrix DoEntrywise(IMatrixView matrix, Func<double, double, double> binaryOperation)
		{
			if (matrix is Matrix dense) return DoEntrywise(dense, binaryOperation);
			else return matrix.DoEntrywise(this, (x, y) => binaryOperation(y, x)); // To avoid accessing zero entries.
		}

		/// <summary>
		/// See <see cref="IEntrywiseOperableView2D{TMatrixIn, TMatrixOut}.DoEntrywise(TMatrixIn, Func{double, double, double})"/>.
		/// </summary>
		public Matrix DoEntrywise(Matrix matrix, Func<double, double, double> binaryOperation)
		{
			Preconditions.CheckSameMatrixDimensions(this, matrix);
			var result = new double[data.Length];
			for (int i = 0; i < data.Length; ++i) result[i] = binaryOperation(this.data[i], matrix.data[i]);
			return new Matrix(result, NumRows, NumColumns);
		}

		/// <summary>
		/// See <see cref="IEntrywiseOperable2D{TMatrixIn}.DoEntrywiseIntoThis(TMatrixIn, Func{double, double, double})"/>.
		/// </summary>
		public void DoEntrywiseIntoThis(IMatrixView matrix, Func<double, double, double> binaryOperation)
		{
			if (matrix is Matrix dense) DoEntrywiseIntoThis(dense, binaryOperation);
			else
			{
				Preconditions.CheckSameMatrixDimensions(this, matrix);
				for (int j = 0; j < NumColumns; ++j)
				{
					for (int i = 0; i < NumRows; ++i)
					{
						int index1D = j * NumRows + i;
						this.data[index1D] = binaryOperation(this.data[index1D], matrix[i, j]);
					}
				}
			}
		}

		/// <summary>
		/// See <see cref="IEntrywiseOperable2D{TMatrixIn}.DoEntrywiseIntoThis(TMatrixIn, Func{double, double, double})"/>.
		/// </summary>
		public void DoEntrywiseIntoThis(Matrix matrix, Func<double, double, double> binaryOperation)
		{
			Preconditions.CheckSameMatrixDimensions(this, matrix);
			for (int i = 0; i < data.Length; ++i) this.data[i] = binaryOperation(this.data[i], matrix.data[i]);
		}

		/// <summary>
		/// See <see cref="IEntrywiseOperableView2D{TMatrixIn, TMatrixOut}.DoToAllEntries(Func{double, double})"/>.
		/// </summary>
		IMatrix IEntrywiseOperableView2D<IMatrixView, IMatrix>.DoToAllEntries(Func<double, double> unaryOperation)
			=> DoToAllEntries(unaryOperation);

		/// <summary>
		/// See <see cref="IEntrywiseOperableView2D{TMatrixIn, TMatrixOut}.DoToAllEntries(Func{double, double})"/>.
		/// </summary>
		public Matrix DoToAllEntries(Func<double, double> unaryOperation)
		{
			var result = new double[NumRows * NumColumns];
			for (int i = 0; i < NumRows * NumColumns; ++i) result[i] = unaryOperation(data[i]);
			return new Matrix(result, NumRows, NumColumns);
		}

		/// <summary>
		/// See <see cref="IEntrywiseOperable2D{TMatrixIn}.DoToAllEntriesIntoThis(Func{double, double})"/>.
		/// </summary>
		public void DoToAllEntriesIntoThis(Func<double, double> unaryOperation)
		{
			for (int i = 0; i < NumRows * NumColumns; ++i) data[i] = unaryOperation(data[i]);
		}

		/// <summary>
		/// See <see cref="IIndexable2D.Equals(IIndexable2D, double)"/>.
		/// </summary>
		public bool Equals(IIndexable2D other, double tolerance = 1e-13)
		{
			if (other is Matrix dense)
			{
				//Check each dimension, rather than the lengths of the internal buffers
				if (!Preconditions.AreSameMatrixDimensions(this, dense)) return false;
				double[] otherData = dense.data;
				var comparer = new ValueComparer(tolerance);
				for (int i = 0; i < this.data.Length; ++i)
				{
					if (!comparer.AreEqual(this.data[i], otherData[i])) return false;
				}
				return true;
			}
			else return other.Equals(this, tolerance); // To avoid accessing zero entries
		}

		/// <summary>
		/// Calculates the Cholesky factorization of a symmetric positive definite matrix with n = <see cref="NumRows"/> = 
		/// <see cref="NumColumns"/>, such that A = L^T * L. L is a lower triangular n-by-n matrix. This only works if the matrix
		/// is symmetric positive definite. Requires extra available memory n^2 entries. 
		/// </summary>
		/// <param name="inPlace">
		/// False, to copy the internal array before factorization. True, to overwrite it with the factorized data, thus saving 
		/// memory and time. However, that will make this object unusable, so you MUST NOT call any other members afterwards.
		/// </param>
		/// <exception cref="NonMatchingDimensionsException">Thrown if the matrix is not square.</exception>
		/// <exception cref="IndefiniteMatrixException">Thrown if the matrix is not symmetric positive definite.</exception>
		/// <exception cref="LapackException">Thrown if the call to LAPACK fails due to invalid input.</exception>
		public CholeskyFull FactorCholesky(bool inPlace = false)
		{
			Preconditions.CheckSquare(this);
			if (inPlace)
			{
				var factor = CholeskyFull.Factorize(NumColumns, data);
				// Set the internal array to null to force NullReferenceException if it is accessed again.
				// TODO: perhaps there is a better way to handle this.
				data = null;
				isOverwritten = true;
				return factor;
			}
			else return CholeskyFull.Factorize(NumColumns, CopyInternalData());
		}

		/// <summary>
		/// Calculates the LQ factorization of a matrix with m = <see cref="NumRows"/> &lt;= <see cref="NumColumns"/> = n, such 
		/// that A = L * Q. Q is an orthogonal n-by-n matrix and L is a lower trapezoidal m-by-n matrix. Requires extra available  
		/// memory form * n + min(m, n) entries.
		/// </summary>
		/// <param name="inPlace">
		/// False, to copy the internal array before factorization. True, to overwrite it with the factorized data, thus saving 
		/// memory and time. However, that will make this object unusable, so you MUST NOT call any other members afterwards.
		/// </param>
		/// <exception cref="LapackException">Thrown if the call to LAPACK fails due to invalid input.</exception>
		public LQFactorization FactorLQ(bool inPlace = false)
		{
			if (inPlace)
			{
				var factor = LQFactorization.Factorize(NumRows, NumColumns, data);
				// Set the internal array to null to force NullReferenceException if it is accessed again.
				// TODO: perhaps there is a better way to handle this.
				data = null;
				isOverwritten = true;
				return factor;
			}
			else return LQFactorization.Factorize(NumRows, NumColumns, CopyInternalData());
		}

		/// <summary>
		/// Calculates the LUP factorization of a square matrix with n = <see cref="NumRows"/> = <see cref="NumColumns"/>, such 
		/// that A = P * L * U. L is a lower triangular n-by-n matrix. U is an upper triangular n-by-n matrix. P is an n-by-n
		/// permutation matrix. Requires extra available memory n^2 + n entries. 
		/// </summary>
		/// <param name="inPlace">
		/// False, to copy the internal array before factorization. True, to overwrite it with the factorized data, thus saving 
		/// memory and time. However, that will make this object unusable, so you MUST NOT call any other members afterwards.
		/// </param>
		/// <exception cref="NonMatchingDimensionsException">Thrown if the matrix is not square.</exception>
		/// <exception cref="LapackException">Thrown if the call to LAPACK fails due to invalid input.</exception>
		public LUFactorization FactorLU(bool inPlace = false)
		{
			Preconditions.CheckSquare(this);
			if (inPlace)
			{
				var factor = LUFactorization.Factorize(NumColumns, data);
				// Set the internal array to null to force NullReferenceException if it is accessed again.
				// TODO: perhaps there is a better way to handle this.
				data = null;
				isOverwritten = true;
				return factor;
			}
			else return LUFactorization.Factorize(NumColumns, CopyInternalData());
		}

		/// <summary>
		/// Calculates the QR factorization of a matrix with m = <see cref="NumRows"/> &gt;= <see cref="NumColumns"/> = n, such 
		/// that A = Q * R. Q is an orthogonal m-by-m matrix and R is an upper trapezoidal m-by-n matrix. Requires extra 
		/// available memory for m * n + min(m, n) entries. 
		/// </summary>
		/// <param name="inPlace">
		/// False, to copy the internal array before factorization. True, to overwrite it with the factorized data, thus saving 
		/// memory and time. However, that will make this object unusable, so you MUST NOT call any other members afterwards.
		/// </param>
		/// <exception cref="LapackException">Thrown if the call to LAPACK fails due to invalid input.</exception>
		public QRFactorization FactorQR(bool inPlace = false)
		{
			if (inPlace)
			{
				var factor = QRFactorization.Factorize(NumRows, NumColumns, data);
				// Set the internal array to null to force NullReferenceException if it is accessed again.
				// TODO: perhaps there is a better way to handle this.
				data = null;
				isOverwritten = true;
				return factor;
			}
			else return QRFactorization.Factorize(NumRows, NumColumns, CopyInternalData());
		}

		/// <summary>
		/// See <see cref="ISliceable2D.GetColumn(int)"/>.
		/// </summary>
		public Vector GetColumn(int colIndex)
		{
			if (isOverwritten) throw new MatrixDataOverwrittenException();
			Preconditions.CheckIndexCol(this, colIndex);
			double[] columnVector = new double[NumRows];
			Array.Copy(data, colIndex * NumRows, columnVector, 0, NumRows);
			return Vector.CreateFromArray(columnVector, false);
		}

		/// <summary>
		/// Returns a <see cref="Vector"/> with the entries of the matrix's main diagonal. The matrix must be square.
		/// </summary>
		/// <exception cref="NonMatchingDimensionsException">Thrown if the matrix is not square.</exception>
		public Vector GetDiagonal() => Vector.CreateFromArray(GetDiagonalAsArray(), false);

		/// <summary>
		/// Returns an array with the entries of the matrix's main diagonal. The matrix must be square.
		/// </summary>
		/// <exception cref="NonMatchingDimensionsException">Thrown if the matrix is not square.</exception>
		public double[] GetDiagonalAsArray()
		{
			if (isOverwritten) throw new MatrixDataOverwrittenException();
			return ArrayColMajor.DiagonalGet(NumRows, NumColumns, data);
		}

		/// <summary>
		/// See <see cref="ISliceable2D.GetRow(int)"/>.
		/// </summary>
		public Vector GetRow(int rowIndex)
		{
			if (isOverwritten) throw new MatrixDataOverwrittenException();
			Preconditions.CheckIndexRow(this, rowIndex);
			double[] rowVector = new double[NumColumns];
			for (int j = 0; j < NumColumns; ++j) rowVector[j] = data[j * NumRows + rowIndex];
			return Vector.CreateFromArray(rowVector, false);
		}

		/// <summary>
		/// See <see cref="ISliceable2D.GetSubmatrix(int[], int[])"/>.
		/// </summary>
		IMatrix ISliceable2D.GetSubmatrix(int[] rowIndices, int[] colIndices) => GetSubmatrix(rowIndices, colIndices);

		public Matrix GetSubmatrix(int[] rowIndices, int[] colIndices)
		{
			if (isOverwritten) throw new MatrixDataOverwrittenException();
			double[] submatrix = new double[colIndices.Length * rowIndices.Length];
			int idxCounter = -1;
			foreach (var j in colIndices)
			{
				foreach (var i in rowIndices) submatrix[++idxCounter] = data[j * NumRows + i];
			}
			return new Matrix(submatrix, rowIndices.Length, colIndices.Length);
		}

		/// <summary>
		/// See <see cref="ISliceable2D.GetSubmatrix(int, int, int, int)"/>.
		/// </summary>
		IMatrix ISliceable2D.GetSubmatrix(int rowStartInclusive, int rowEndExclusive, int colStartInclusive, int colEndExclusive)
			=> GetSubmatrix(rowStartInclusive, rowEndExclusive, colStartInclusive, colEndExclusive);

		public Matrix GetSubmatrix(int rowStartInclusive, int rowEndExclusive, int colStartInclusive, int colEndExclusive)
		{
			if (isOverwritten) throw new MatrixDataOverwrittenException();
			int numNewRows = rowEndExclusive - rowStartInclusive;
			int numNewCols = colEndExclusive - colStartInclusive;
			double[] submatrix = new double[numNewCols * numNewRows];
			int idxCounter = -1;
			for (int j = colStartInclusive; j < colEndExclusive; ++j)
			{
				for (int i = rowStartInclusive; i < rowEndExclusive; ++i)
				{
					submatrix[++idxCounter] = data[j * NumRows + i];
				}
			}
			return new Matrix(submatrix, numNewRows, numNewCols);
		}

		/// <summary>
		/// Calculates the inverse matrix and returns it in a new <see cref="Matrix"/> instance. This only works if this 
		/// <see cref="Matrix"/> is square and invertible. If the determinant matrix is also needed, use 
		/// <see cref="InvertAndDeterminant"/> instead.
		/// </summary>
		/// <exception cref="NonMatchingDimensionsException">Thrown if the matrix is not square.</exception>
		/// <exception cref="SingularMatrixException">Thrown if the matrix is not invertible.</exception>
		/// <exception cref="LapackException">Thrown if the call to LAPACK fails due to invalid input.</exception>
		public Matrix Invert() => Invert(AnalyticFormulas.determinantTolerance);

		/// <summary>
		/// Calculates the inverse matrix and returns it in a new <see cref="Matrix"/> instance. This only works if this 
		/// <see cref="Matrix"/> is square and invertible. If the determinant matrix is also needed, use 
		/// <see cref="InvertAndDeterminant"/> instead.
		/// </summary>
		/// <param name="tolerance">
		/// Matrix determinant value under which the matrix is considered singular (valid ONLY for 2x2 and 3x3 matrices)
		/// </param>
		/// <exception cref="NonMatchingDimensionsException">Thrown if the matrix is not square.</exception>
		/// <exception cref="SingularMatrixException">Thrown if the matrix is not invertible.</exception>
		/// <exception cref="LapackException">Thrown if the call to LAPACK fails due to invalid input.</exception>
		public Matrix Invert(double tolerance)
		{
			if ((NumRows == 2) && (NumColumns == 2))
			{
				(double[] inverse, double det) = AnalyticFormulas.Matrix2x2ColMajorInvert(data, tolerance);
				return new Matrix(inverse, 2, 2);
			}
			else if ((NumRows == 3) && (NumColumns == 3))
			{
				(double[] inverse, double det) = AnalyticFormulas.Matrix3x3ColMajorInvert(data, tolerance);
				return new Matrix(inverse, 3, 3);
			}
			else
			{
				return FactorLU(false).Invert(true);
			}
		}

		/// <summary>
		/// Calculates the inverse matrix and writes it over the entries of this object, in order to conserve memory 
		/// and possibly time.
		/// </summary>
		/// <exception cref="NonMatchingDimensionsException">Thrown if the matrix is not square.</exception>
		/// <exception cref="SingularMatrixException">Thrown if the matrix is not invertible.</exception>
		/// <exception cref="LapackException">Thrown if the call to LAPACK fails due to invalid input.</exception>
		public void InvertInPlace() => InvertInPlace(AnalyticFormulas.determinantTolerance);

		/// <summary>
		/// Calculates the inverse matrix and writes it over the entries of this object, in order to conserve memory 
		/// and possibly time.
		/// </summary>
		/// <param name="tolerance">
		/// Matrix determinant value under which the matrix is considered singular (valid ONLY for 2x2 and 3x3 matrices)
		/// </param>
		/// <exception cref="NonMatchingDimensionsException">Thrown if the matrix is not square.</exception>
		/// <exception cref="SingularMatrixException">Thrown if the matrix is not invertible.</exception>
		/// <exception cref="LapackException">Thrown if the call to LAPACK fails due to invalid input.</exception>
		public void InvertInPlace(double tolerance) //TODO: This needs redesign. A new object should be returned and this should be disabled.
		{
			//TODO: implement efficient 2x2 and 3x3 inplace operations to avoid copying. 
			if ((NumRows == 2) && (NumColumns == 2))
			{
				(double[] inverse, double det) = AnalyticFormulas.Matrix2x2ColMajorInvert(data, tolerance);
				Array.Copy(inverse, data, 4);
			}
			else if ((NumRows == 3) && (NumColumns == 3))
			{
				(double[] inverse, double det) = AnalyticFormulas.Matrix3x3ColMajorInvert(data, tolerance);
				Array.Copy(inverse, data, 9);
			}
			else
			{
				// The next will update the entries of this matrix, but we do not need the intermediate objects
				LUFactorization.Factorize(NumColumns, data).Invert(true); 
			}
		}

		/// <summary>
		/// Calculates the determinant and the inverse matrix and returns the latter in a new <see cref="Matrix"/> instance. 
		/// This only works if this <see cref="Matrix"/> is square and invertible.
		/// </summary>
		/// <exception cref="NonMatchingDimensionsException">Thrown if the matrix is not square.</exception>
		/// <exception cref="SingularMatrixException">Thrown if the matrix is not invertible.</exception>
		/// <exception cref="LapackException">Thrown if the call to LAPACK fails due to invalid input.</exception>
		public (Matrix inverse, double determinant) InvertAndDeterminant() => InvertAndDeterminant(AnalyticFormulas.determinantTolerance);

		/// <summary>
		/// Calculates the determinant and the inverse matrix and returns the latter in a new <see cref="Matrix"/> instance. 
		/// This only works if this <see cref="Matrix"/> is square and invertible.
		/// </summary>
		/// <param name="tolerance">
		/// Matrix determinant value under which the matrix is considered singular (valid ONLY for 2x2 and 3x3 matrices)
		/// </param>
		/// <exception cref="NonMatchingDimensionsException">Thrown if the matrix is not square.</exception>
		/// <exception cref="SingularMatrixException">Thrown if the matrix is not invertible.</exception>
		/// <exception cref="LapackException">Thrown if the call to LAPACK fails due to invalid input.</exception>
		public (Matrix inverse, double determinant) InvertAndDeterminant(double tolerance)
		{
			if ((NumRows == 2) && (NumColumns == 2))
			{
				(double[] inverse, double det) = AnalyticFormulas.Matrix2x2ColMajorInvert(data, tolerance);
				return (new Matrix(inverse, 2, 2), det);
			}
			else if ((NumRows == 3) && (NumColumns == 3))
			{
				(double[] inverse, double det) = AnalyticFormulas.Matrix3x3ColMajorInvert(data, tolerance);
				return (new Matrix(inverse, 3, 3), det);
			}
			else
			{
				LUFactorization factor = FactorLU(false);
				double det = factor.CalcDeterminant(); // Call this before factor.Invert(), else the factor will be overwritten.
				Matrix inverse = factor.Invert(true); 
				return (inverse, det);
			}
		}

		/// <summary>
		/// Returns true if: this[i, j] &lt;= <paramref name="tolerance"/>, for 0 &lt;= i &lt; <see cref="NumRows"/>, 
		/// 0 &lt;= j &lt; <see cref="NumColumns"/>. Otherwise false is returned.
		/// </summary>
		/// <param name="tolerance">The tolerance under which a matrix entry is considered to be 0. It can be set to 0, to check 
		///     if the entries are exactly 0.</param>
		public bool IsZero(double tolerance) => DenseStrategies.IsZero(data, tolerance);

		/// <summary>
		/// See <see cref="IMatrixView.LinearCombination(double, IMatrixView, double)"/>.
		/// </summary>
		public IMatrix LinearCombination(double thisCoefficient, IMatrixView otherMatrix, double otherCoefficient)
		{
			if (otherMatrix is Matrix dense) return LinearCombination(thisCoefficient, dense, otherCoefficient);
			else return otherMatrix.LinearCombination(otherCoefficient, this, thisCoefficient); // To avoid accessing zero entries
		}

		/// <summary>
		/// Performs the following operation for 0 &lt;= i &lt; <see cref="NumRows"/>, 0 &lt;= j &lt; <see cref="NumColumns"/>:
		/// result[i, j] = <paramref name="thisCoefficient"/> * this[i, j] 
		///     + <paramref name="otherCoefficient"/> * <paramref name="otherMatrix"/>[i, j]. 
		/// The resulting matrix is written to a new <see cref="Matrix"/> and then returned.
		/// </summary>
		/// <param name="thisCoefficient">A scalar that multiplies each entry of this <see cref="Matrix"/>.</param>
		/// <param name="otherMatrix">A matrix with the same <see cref="NumRows"/> and <see cref="NumColumns"/> as this 
		///     <see cref="Matrix"/> instance.</param>
		/// <param name="otherCoefficient">A scalar that multiplies each entry of <paramref name="otherMatrix"/>.</param>
		/// <exception cref="NonMatchingDimensionsException">Thrown if <paramref name="otherMatrix"/> has different 
		///     <see cref="NumRows"/> or <see cref="NumColumns"/> than this instance.</exception>
		public Matrix LinearCombination(double thisCoefficient, Matrix otherMatrix, double otherCoefficient)
		{
			Preconditions.CheckSameMatrixDimensions(this, otherMatrix);
			//TODO: Perhaps this should be done using mkl_malloc and BLAS copy. 
			double[] result = new double[data.Length];
			if (thisCoefficient == 1.0)
			{
				Array.Copy(this.data, result, data.Length);
				Blas.Daxpy(data.Length, otherCoefficient, otherMatrix.data, 0, 1, result, 0, 1);
			}
			else if (otherCoefficient == 1.0)
			{
				Array.Copy(otherMatrix.data, result, data.Length);
				Blas.Daxpy(data.Length, thisCoefficient, this.data, 0, 1, result, 0, 1);
			}
			else
			{
				Array.Copy(this.data, result, data.Length);
				BlasExtensions.Daxpby(data.Length, otherCoefficient, otherMatrix.data, 0, 1, thisCoefficient, result, 0, 1);
			}
			return new Matrix(result, NumRows, NumColumns);
		}

		/// <summary>
		/// See <see cref="IMatrix.LinearCombinationIntoThis(double, IMatrixView, double)"/>.
		/// </summary>
		public void LinearCombinationIntoThis(double thisCoefficient, IMatrixView otherMatrix, double otherCoefficient)
		{
			if (otherMatrix is Matrix dense) LinearCombinationIntoThis(thisCoefficient, dense, otherCoefficient);
			else
			{
				Preconditions.CheckSameMatrixDimensions(this, otherMatrix);
				for (int j = 0; j < NumColumns; ++j)
				{
					for (int i = 0; i < NumRows; ++i)
					{
						int index1D = j * NumRows + i;
						this.data[index1D] = thisCoefficient * this.data[index1D] + otherCoefficient * otherMatrix[i, j];
					}
				}
			}
		}

		/// <summary>
		/// Performs the following operation for 0 &lt;= i &lt; <see cref="NumRows"/>, 0 &lt;= j &lt; <see cref="NumColumns"/>:
		/// this[i, j] = <paramref name="thisCoefficient"/> * this[i, j] 
		///     + <paramref name="otherCoefficient"/> * <paramref name="otherMatrix"/>[i, j]. 
		/// The resulting matrix overwrites the entries of this <see cref="Matrix"/> instance.
		/// </summary>
		/// <param name="thisCoefficient">A scalar that multiplies each entry of this <see cref="Matrix"/>.</param>
		/// <param name="otherMatrix">A matrix with the same <see cref="NumRows"/> and <see cref="NumColumns"/> as this 
		///     <see cref="Matrix"/> instance.</param>
		/// <param name="otherCoefficient">A scalar that multiplies each entry of <paramref name="otherMatrix"/>.</param>
		/// <exception cref="NonMatchingDimensionsException">Thrown if <paramref name="otherMatrix"/> has different 
		///     <see cref="NumRows"/> or <see cref="NumColumns"/> than this instance.</exception>
		public void LinearCombinationIntoThis(double thisCoefficient, Matrix otherMatrix, double otherCoefficient)
		{
			Preconditions.CheckSameMatrixDimensions(this, otherMatrix);
			if (thisCoefficient == 1.0)
			{
				Blas.Daxpy(data.Length, otherCoefficient, otherMatrix.data, 0, 1, this.data, 0, 1);
			}
			else
			{
				BlasExtensions.Daxpby(data.Length, otherCoefficient, otherMatrix.data, 0, 1, thisCoefficient, this.data, 0, 1);
			}
		}

		/// <summary>
		/// See <see cref="IMatrixView.MultiplyLeft(IMatrixView, bool, bool)"/>.
		/// </summary>
		public Matrix MultiplyLeft(IMatrixView other, bool transposeThis = false, bool transposeOther = false)
		{
			return other.MultiplyRight(this, transposeOther, transposeThis);
		}

		/// <summary>
		/// See <see cref="IMatrixView.MultiplyRight(IMatrixView, bool, bool)"/>.
		/// </summary>
		public Matrix MultiplyRight(IMatrixView other, bool transposeThis = false, bool transposeOther = false)
		{
			if (other is Matrix dense) return MultiplyRight(dense, transposeThis, transposeOther);
			else return other.MultiplyLeft(this, transposeOther, transposeThis);
		}

		/// <summary>
		/// Performs the matrix-matrix multiplication: oper(this) * oper(<paramref name="other"/>).
		/// </summary>
		/// <param name="other">A matrix such that the <see cref="NumRows"/> of oper(<paramref name="other"/>) 
		///     are equal to the <see cref="NumColumns"/> of oper(this).</param>
		/// <param name="transposeThis">If true, oper(this) = transpose(this). Otherwise oper(this) = this.</param>
		/// <param name="transposeOther">If true, oper(<paramref name="other"/>) = transpose(<paramref name="other"/>). 
		///     Otherwise oper(<paramref name="other"/>) = <paramref name="other"/>.</param>
		/// <exception cref="NonMatchingDimensionsException">Thrown if oper(<paramref name="otherMatrix"/>) has 
		///     different <see cref="NumRows"/> than the <see cref="NumColumns"/> of oper(this).</exception>
		public Matrix MultiplyRight(Matrix other, bool transposeThis = false, bool transposeOther = false)
		{
			int leftRows, leftCols, rightRows, rightCols;
			TransposeMatrix transposeLeft, transposeRight;
			if (transposeThis)
			{
				transposeLeft = TransposeMatrix.Transpose;
				leftRows = this.NumColumns;
				leftCols = this.NumRows;
			}
			else
			{
				transposeLeft = TransposeMatrix.NoTranspose;
				leftRows = this.NumRows;
				leftCols = this.NumColumns;
			}
			if (transposeOther)
			{
				transposeRight = TransposeMatrix.Transpose;
				rightRows = other.NumColumns;
				rightCols = other.NumRows;
			}
			else
			{
				transposeRight = TransposeMatrix.NoTranspose;
				rightRows = other.NumRows;
				rightCols = other.NumColumns;
			}

			Preconditions.CheckMultiplicationDimensions(leftCols, rightRows);
			double[] result = new double[leftRows * rightCols];
			Blas.Dgemm(transposeLeft, transposeRight, leftRows, rightCols, leftCols,
				1.0, this.data, 0, this.NumRows, other.data, 0, other.NumRows,
				1.0, result, 0, leftRows);
			return new Matrix(result, leftRows, rightCols);
		}

		/// <summary>
		/// See <see cref="IMatrixView.Multiply(IVectorView, bool)"/>.
		/// </summary>
		public IVector Multiply(IVectorView vector, bool transposeThis = false)
		{
			if (vector is Vector dense) return Multiply(dense, transposeThis);
			else throw new NotImplementedException();
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

		/// <summary>
		/// See <see cref="IMatrixView.MultiplyIntoResult(IVectorView, IVector, bool)"/>.
		/// </summary>
		public void MultiplyIntoResult(IVectorView lhsVector, IVector rhsVector, bool transposeThis = false)
		{
			if ((lhsVector is Vector lhsDense) && (rhsVector is Vector rhsDense))
			{
				MultiplyIntoResult(lhsDense, rhsDense, transposeThis);
			}
			else throw new NotImplementedException();
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
		{   //TODO: this is NOT a specialization of the version with offsets. It is defined only if the vectors have exactly the matching lengths.
			int leftRows, leftCols;
			TransposeMatrix transpose;
			if (transposeThis)
			{
				transpose = TransposeMatrix.Transpose;
				leftRows = this.NumColumns;
				leftCols = this.NumRows;
			}
			else
			{
				transpose = TransposeMatrix.NoTranspose;
				leftRows = this.NumRows;
				leftCols = this.NumColumns;
			}

			Preconditions.CheckMultiplicationDimensions(leftCols, lhsVector.Length);
			Preconditions.CheckSystemSolutionDimensions(leftRows, rhsVector.Length);
			Blas.Dgemv(transpose, NumRows, NumColumns,
				1.0, this.data, 0, NumRows, lhsVector.RawData, 0, 1,
				0.0, rhsVector.RawData, 0, 1);
		}

		/// <summary>
		/// Performs the matrix-vector multiplication y = alpha * oper(A) * x + beta * y, where:
		/// alpha = <paramref name="lhsScale"/>, x = <paramref name="lhsVector"/>, beta = <paramref name="rhsScale"/>,
		/// y = <paramref name="rhsScale"/>, oper(A) = this (if <paramref name="transposeThis"/> == false) or transpose(this)
		/// (if <paramref name="transposeThis"/> == true).
		/// The input vectors can be longer (taking into account the offsets) than the corresponding dimensions of this matrix. 
		/// The resulting vector will overwrite <paramref name="lhsVector"/> starting from <paramref name="rhsOffset"/>.
		/// </summary>
		/// <param name="lhsVector">
		/// The vector x that will be multiplied by this matrix. Constraints: 
		/// <paramref name="lhsOffset"/> + <paramref name="lhsVector"/>.<see cref="IIndexable1D.Length"/> 
		/// &lt;= oper(this).<see cref="IIndexable2D.NumColumns"/>.
		/// </param>
		/// <param name="lhsOffset">The index into <paramref name="lhsVector"/> from which to start the operations.</param>
		/// <param name="lhsScale">The scalar alpha that will multiply <paramref name="lhsVector"/>.</param>
		/// <param name="rhsVector">
		/// The vector y that will be overwritten by the result of the operation. Constraints: 
		/// <paramref name="rhsOffset"/> + <paramref name="rhsVector"/>.<see cref="IIndexable1D.Length"/> 
		/// &lt;= oper(this).<see cref="IIndexable2D.NumRows"/>.
		/// </param>
		/// <param name="rhsOffset">The index into <paramref name="rhsVector"/> from which to start the operations.</param>
		/// <param name="rhsScale">The scalar beta that will multiply <paramref name="rhsVector"/>.</param>
		/// <param name="transposeThis">If true, oper(this) = transpose(this). Otherwise oper(this) = this.</param>
		/// <exception cref="NonMatchingDimensionsException">
		/// Thrown if the <see cref="IIndexable1D.Length"/> of <paramref name="lhsVector"/> or <paramref name="rhsVector"/> 
		/// violate the described contraints.
		/// </exception>
		/// <exception cref="PatternModifiedException">
		/// Thrown if the storage format of <paramref name="rhsVector"/> does not support overwritting the entries that this 
		/// method will try to.
		/// </exception>
		public void MultiplySubvectorIntoResult(Vector lhsVector, int lhsOffset, double lhsScale, Vector rhsVector, int rhsOffset,
			double rhsScale, bool transposeThis = false)
		{
			Preconditions.CheckMultiplicationDimensions(this, lhsVector, lhsOffset, rhsVector, rhsOffset, transposeThis);
			if (transposeThis)
			{
				Blas.Dgemv(TransposeMatrix.Transpose, NumColumns, NumRows, lhsScale, this.data, 0, NumRows,
					lhsVector.RawData, lhsOffset, 1, rhsScale, rhsVector.RawData, rhsOffset, 1);
			}
			else
			{
				Blas.Dgemv(TransposeMatrix.NoTranspose, NumRows, NumColumns, lhsScale, this.data, 0, NumRows,
					lhsVector.RawData, lhsOffset, 1, rhsScale, rhsVector.RawData, rhsOffset, 1);
			}
		}

		/// <summary>
		/// See <see cref="IReducible.Reduce(double, ProcessEntry, ProcessZeros, Reduction.Finalize)"/>.
		/// </summary>
		public double Reduce(double identityValue, ProcessEntry processEntry, ProcessZeros processZeros, Finalize finalize)
		{
			double aggregator = identityValue;
			for (int i = 0; i < data.Length; ++i) aggregator = processEntry(data[i], aggregator);
			// no zeros implied
			return finalize(aggregator);
		}

		/// <summary>
		/// Creates a new <see cref="Matrix"/> that contains the columns of this <see cref="Matrix"/> with a different order,
		/// which is specified by the provided <paramref name="permutation"/> and <paramref name="oldToNew"/>.
		/// </summary>
		/// <param name="permutation">
		/// An array that contains the row/column indices of this <see cref="Matrix"/> in a different order.
		/// </param>
		/// <param name="oldToNew">
		/// If true, reordered[i, <paramref name="permutation"/>[j]] =  original[i, j]. 
		/// If false, reordered[i, j] = original[i, <paramref name="permutation"/>[j]].
		/// </param>
		/// <exception cref="NonMatchingDimensionsException">
		/// If <paramref name="permutation"/>.Length is different than the number of columns of this matrix.
		/// </exception>
		public Matrix ReorderColumns(IReadOnlyList<int> permutation, bool oldToNew)
		{
			if (permutation.Count != NumColumns) throw new NonMatchingDimensionsException(
				$"This matrix has {NumColumns} columns, while the permutation vector has {permutation.Count} entries.");
			if (oldToNew)
			{
				return new Matrix(ArrayColMajor.ReorderColumnsOldToNew(NumRows, NumColumns, data, permutation), NumRows, NumRows);
			}
			else
			{
				return new Matrix(ArrayColMajor.ReorderColumnsNewToOld(NumRows, NumColumns, data, permutation), NumRows, NumRows);
			}
		}

		/// <summary>
		/// Creates a new <see cref="Matrix"/> that contains the entries of this <see cref="Matrix"/> with a different order,
		/// which is specified by the provided <paramref name="permutation"/> and <paramref name="oldToNew"/>.
		/// </summary>
		/// <param name="permutation">
		/// An array that contains the row/column indices of this <see cref="Matrix"/> in a different order.
		/// </param>
		/// <param name="oldToNew">
		/// If true, reordered[<paramref name="permutation"/>[i], <paramref name="permutation"/>[j]] =  original[i, j]. 
		/// If false, reordered[i, j] = original[<paramref name="permutation"/>[i], <paramref name="permutation"/>[j]].
		/// </param>
		/// /// <exception cref="NonMatchingDimensionsException">
		/// If this matrix is not square or its number of rows/columns if different than <paramref name="permutation"/>.Length.
		/// </exception>
		public Matrix Reorder(IReadOnlyList<int> permutation, bool oldToNew)
		{
			Preconditions.CheckSquare(this);
			if (permutation.Count != NumRows) throw new NonMatchingDimensionsException(
				$"This matrix has order = {NumRows}, while the permutation vector has {permutation.Count} entries.");
			if (oldToNew) return new Matrix(ArrayColMajor.ReorderOldToNew(NumRows, data, permutation), NumRows, NumRows);
			else return new Matrix(ArrayColMajor.ReorderNewToOld(NumRows, data, permutation), NumRows, NumRows);
		}

		/// <summary>
		/// See <see cref="IMatrixView.Scale(double)"/>.
		/// </summary>
		IMatrix IMatrixView.Scale(double scalar) => Scale(scalar);

		/// <summary>
		/// Performs the following operation for 0 &lt;= i &lt; <see cref="NumRows"/>, 0 &lt;= j &lt; <see cref="NumColumns"/>:
		/// result[i, j] = <paramref name="scalar"/> * this[i, j].
		/// The resulting matrix is written to a new <see cref="Matrix"/> and then returned.
		/// </summary>
		/// <param name="scalar">A scalar that multiplies each entry of this matrix.</param>
		public Matrix Scale(double scalar)
		{
			//TODO: Perhaps this should be done using mkl_malloc and BLAS copy. 
			double[] result = new double[data.Length];
			Array.Copy(data, result, data.Length);
			Blas.Dscal(data.Length, scalar, result, 0, 1);
			return new Matrix(result, NumRows, NumColumns);
		}

		/// <summary>
		/// See <see cref="IMatrix.ScaleIntoThis(double)"/>.
		/// </summary>
		public void ScaleIntoThis(double scalar) => Blas.Dscal(data.Length, scalar, data, 0, 1);

		/// <summary>
		/// Sets all entries of this matrix to be equal to <paramref name="value"/>.
		/// </summary>
		/// <param name="value">The value that all entries of the this matrix will be equal to.</param>
		public void SetAll(double value)
		{
			for (int i = 0; i < data.Length; ++i) data[i] = value;
		}

		/// <summary>
		/// Sets some consecutive entries of the column with index = <paramref name="colIdx"/> to be equal to 
		/// <paramref name="colValues"/>, starting from the entry with row index = <paramref name="rowStart"/>.
		/// </summary>
		/// <param name="colIdx">The index of the column to set. Constraints: 
		///     0 &lt;= <paramref name="colIdx"/> &lt; <see cref="NumColumns"/>.</param>
		/// <param name="rowStart">The first entry of column <paramref name="colIdx"/> to be modified. Constraints: 
		///     1) 0 &lt;= <paramref name="rowStart"/> &lt; <see cref="NumRows"/>, 
		///     2) <paramref name="rowStart"/> + <paramref name="colValues"/>.<see cref="IIndexable1D.Length"/> &lt;= 
		///        <see cref="NumRows"/>.</param>
		/// <param name="colValues">The new values of the column entries. Constraints: <paramref name="rowStart"/>
		///     + <paramref name="colValues"/>.<see cref="IIndexable1D.Length"/> &lt;= <see cref="NumRows"/>.</param>
		/// <exception cref="IndexOutOfRangeException">Thrown if <paramref name="colIdx"/> or <paramref name="rowStart"/> 
		///     violate the described constraints.</exception>
		/// <exception cref="NonMatchingDimensionsException">Thrown if <paramref name="rowStart"/>
		///     + <paramref name="colValues"/>.<see cref="IIndexable1D.Length"/> &gt; <see cref="NumRows"/>.</exception>
		public void SetSubcolumn(int colIdx, Vector colValues, int rowStart = 0)
		{
			Preconditions.CheckIndexCol(this, colIdx);
			if (rowStart + colValues.Length > this.NumRows) throw new NonMatchingDimensionsException(
				"The entries to set exceed this matrix's number of rows");
			ArrayColMajor.SetCol(NumRows, NumColumns, data, colIdx, rowStart, colValues.RawData);
		}

		/// <summary>
		/// See <see cref="IMatrix.SetEntryRespectingPattern(int, int, double)"/>.
		/// </summary>
		public void SetEntryRespectingPattern(int rowIdx, int colIdx, double value)
		{
			data[colIdx * NumRows + rowIdx] = value;
		}

		/// <summary>
		/// Sets some consecutive entries of this matrix to be equal to the entries of <paramref name="submatrix"/>, starting 
		/// from the entry at (<paramref name="rowStart"/>, <paramref name="colStart"/>).
		/// </summary>
		/// <param name="rowStart">The index of the first row to be modified. Constraints: 
		///     1) 0 &lt;= <paramref name="rowStart"/> &lt; <see cref="NumRows"/>, 
		///     2) <paramref name="rowStart"/> + <paramref name="submatrix"/>.<see cref="NumRows"/> &lt;= 
		///        this.<see cref="NumRows"/>.</param>
		/// <param name="colStart">The index of the first column to be modified. Constraints: 
		///     1) 0 &lt;= <paramref name="colStart"/> &lt; <see cref="NumColumns"/>, 
		///     2) <paramref name="colStart"/> + <paramref name="submatrix"/>.<see cref="NumColumns"/> &lt;= 
		///        this.<see cref="NumColumns"/>.</param>
		/// <param name="submatrix">The new values of this matrix's entries to be modified. Constraints:
		///     1) <paramref name="rowStart"/> + <paramref name="submatrix"/>.<see cref="NumRows"/> &lt;= 
		///        this.<see cref="NumRows"/>,
		///     2) <paramref name="colStart"/> + <paramref name="submatrix"/>.<see cref="NumColumns"/> &lt;= 
		///        this.<see cref="NumColumns"/>.</param>
		/// <exception cref="IndexOutOfRangeException">Thrown if <paramref name="rowStart"/> or <paramref name="colStart"/> 
		///     violate the described constraints.</exception>
		/// <exception cref="NonMatchingDimensionsException">Thrown if <paramref name="rowStart"/>, <paramref name="colStart"/>
		///     or <paramref name="submatrix"/> exceed the desrcibed constraints.</exception>
		public void SetSubmatrix(int rowStart, int colStart, Matrix submatrix)
		{
			// TODO: create Preconditions.CheckOverflow1D() and 2D for such setters.
			Preconditions.CheckIndices(this, rowStart, colStart);
			if ((rowStart + submatrix.NumRows > this.NumRows) || (colStart + submatrix.NumColumns > this.NumColumns))
			{
				throw new NonMatchingDimensionsException("The submatrix doesn't fit inside this matrix, at least when starting"
					+ " from the specified entry.");
			}
			ArrayColMajor.SetSubmatrix(this.NumRows, this.NumColumns, this.data, rowStart, colStart, 
				submatrix.NumRows, submatrix.NumColumns, submatrix.data);
		}

		/// <summary>
		/// Sets some consecutive entries of the row with index = <paramref name="rowIdx"/> to be equal to 
		/// <paramref name="rowValues"/>, starting from the entry with column index = <paramref name="colStart"/>.
		/// </summary>
		/// <param name="rowIdx">The index of the row to set. Constraints: 
		///     0 &lt;= <paramref name="rowIdx"/> &lt; <see cref="NumRows"/>.</param>
		/// <param name="colStart">The first entry of row <paramref name="rowIdx"/> to be modified. Constraints: 
		///     1) 0 &lt;= <paramref name="colStart"/> &lt; <see cref="NumColumns"/>, 
		///     2) <paramref name="colStart"/> + <paramref name="rowValues"/>.<see cref="IIndexable1D.Length"/> &lt;= 
		///        <see cref="NumColumns"/>.</param>
		/// <param name="rowValues">The new values of the row entries. Constraints: <paramref name="colStart"/>
		///     + <paramref name="rowValues"/>.<see cref="IIndexable1D.Length"/> &lt;= <see cref="NumColumns"/>.</param>
		/// <exception cref="IndexOutOfRangeException">Thrown if <paramref name="rowIdx"/> or <paramref name="colStart"/> 
		///     violate the described constraints.</exception>
		/// <exception cref="NonMatchingDimensionsException">Thrown if <paramref name="colStart"/>
		///     + <paramref name="rowValues"/>.<see cref="IIndexable1D.Length"/> &gt; <see cref="NumColumns"/>.</exception>
		public void SetSubrow(int rowIdx, Vector rowValues, int colStart = 0)
		{
			Preconditions.CheckIndexRow(this, rowIdx);
			if (colStart + rowValues.Length > this.NumRows) throw new NonMatchingDimensionsException(
				"The entries to set exceed this matrix's number of columns");
			ArrayColMajor.SetRow(NumRows, NumColumns, data, rowIdx, colStart, rowValues.RawData);
		}

		/// <summary>
		/// Calculates the Singular Value Decomposition of a symmetric matrix.
		/// </summary>
		/// <param name="w">Vector with singular values</param>
		/// <param name="v">Matrix with orthonormal vectors as columns</param>
		public void SVD(double[] w, double[,] v)
		{
			DenseStrategies.SVD(this, w, v);
		}

		/// <summary>
		/// See <see cref="IMatrixView.Transpose"/>.
		/// </summary>
		IMatrix IMatrixView.Transpose() => Transpose();

		/// <summary>
		/// Initializes a new <see cref="Matrix"/> instance, that is transpose to this: result[i, j] = this[j, i]. The entries 
		/// will be explicitly copied.
		/// </summary>
		public Matrix Transpose()
		{
			//TODO: The wrapper library does not include LAPACK's CBlas-like extensions yet. Create my own wrapper or 
			// piggyback on another BLAS function.
			double[] transpose = Conversions.ColumnMajorToRowMajor(data, NumRows, NumColumns);
			return new Matrix(transpose, NumColumns, NumRows);
		}

		/// <summary>
		/// Transposes the matrix by modifying the entries of this <see cref="Matrix instance"/>: this[i, j] = this[j, i].
		/// </summary>
		public void TransposeIntoThis()
		{
			throw new NotImplementedException("Use mkl_dimatcopy");
		}

		private double[] CopyInternalData()
		{
			double[] dataCopy = new double[data.Length];
			Array.Copy(data, dataCopy, data.Length);
			return dataCopy;
		}
	}
}
