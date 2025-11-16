//TODO: Also linear combinations with other matrix types may be useful, e.g. Skyline (K) with diagonal (M), but I think 
//      that for global matrices, this should be done through concrete class to use DoEntrywiseIntoThis methods. 
//TODO: Checks like: col - row <= colHeight can be written more efficiently without calculating the height:
//      entryOffset = diagOffsets[col] + col - row, entryOffset <= diagOffsets[col+1]
//TODO: Throw MatrixDataOverwrittenException whenever necessary. Possibly make the values and diagOffsets readonly again.
//TODO: In most algorithms, I can cache in local variables for each column the height, diagOffset and diagOffset + colIdx
//TODO: Most algorithms implemented here should be moved to a class the holds the implementations and called from there.
namespace MGroup.LinearAlgebra.Matrices
{
	using System;
	using System.Collections.Generic;
	using System.Linq;
	using System.Runtime.CompilerServices;

	using MGroup.LinearAlgebra.Commons;
	using MGroup.LinearAlgebra.Exceptions;
	using MGroup.LinearAlgebra.Extensions;
	using MGroup.LinearAlgebra.Implementations;
	using MGroup.LinearAlgebra.Implementations.Managed;
	using MGroup.LinearAlgebra.Matrices.Builders;
	using MGroup.LinearAlgebra.Output.Formatting;
	using MGroup.LinearAlgebra.Reduction;
	using MGroup.LinearAlgebra.Reordering;
	using MGroup.LinearAlgebra.Triangulation;
	using MGroup.LinearAlgebra.Vectors;

	using static MGroup.LinearAlgebra.LibrarySettings;

	/// <summary>
	/// Symmetric sparse matrix stored in Skyline format (3-array version). Only the non-zero entries of the upper triangle are
	/// stored in column major order. The Skyline format is optimized for Cholesky factorizations.
	/// To build a <see cref="SkylineMatrix"/> conveniently, use <see cref="SkylineBuilder"/>.
	/// </summary>
	public class SkylineMatrix : ValuesBackedMatrix<SkylineMatrix>, ISparseMatrix, ISymmetricMatrix
	{
		/// <summary>
		/// Contains the indices into values of the diagonal entries of the matrix. Its length = order + 1, with the last entry
		/// being equal to nnz.
		/// </summary>
		private int[] diagOffsets;

		/// <summary>
		/// An array containing the nonzero entries of the upper triangle of the matrix.
		/// </summary>
		private double[] values;

		// TODO: This was not a good idea after all. Must check it in every method (very easy to forget). It would be better
		//		 to set values to null.
		private bool isOverwritten = false;

		private SkylineMatrix(int order, double[] values, int[] diagOffsets)
		{
			this.values = values;
			this.diagOffsets = diagOffsets;
			this.NumColumns = order;
		}

		public override MatrixSymmetry MatrixSymmetry => MatrixSymmetry.Symmetric;

		public override int NumColumns { get; }

		public override int NumRows => NumColumns;

		/// <summary>
		/// The internal array that stores the indices into <see cref="RawValues"/> of the diagonal entries of the matrix. 
		/// Its length = order + 1, with the last entry being equal to nnz.
		/// It should only be used for passing the raw array to linear algebra libraries.
		/// </summary>
		public int[] RawDiagOffsets => diagOffsets;

		/// <summary>
		/// The internal array that stores the non-zero entries of the matrix's upper triangle in column major order, 
		/// starting from the diagonal and going upwards. Its length is equal to the number of non-zero entries. 
		/// It should only be used for passing the raw array to linear algebra libraries.
		/// </summary>
		public override double[] RawValues => values;

		public override double this[int rowIdx, int colIdx] // TODO: should I add index bound checking?
		{
			get
			{
				ProcessIndices(ref rowIdx, ref colIdx);
				int diagOffset = diagOffsets[colIdx];
				int maxColumnHeight = diagOffsets[colIdx + 1] - diagOffset - 1; // excluding diagonal
				int entryHeight = colIdx - rowIdx; // excluding diagonal
				if (entryHeight > maxColumnHeight) return 0.0; // outside stored non zero pattern
				else return values[diagOffset + entryHeight];
			}

			set
			{
				ProcessIndices(ref rowIdx, ref colIdx);
				int diagOffset = diagOffsets[colIdx];
				int maxColumnHeight = diagOffsets[colIdx + 1] - diagOffset - 1; // excluding diagonal
				int entryHeight = colIdx - rowIdx; // excluding diagonal
				if (entryHeight > maxColumnHeight)
				{
					throw new SparsityPatternModifiedException($"In column {colIdx} only rows [{maxColumnHeight}, {colIdx}]"
						+ $" can be changed, but you are trying to set entry ({rowIdx}, {colIdx})");
				}
				else values[diagOffset + entryHeight] = value;
			}
		}

		/// <summary>
		/// Initializes a new <see cref="SkylineMatrix"/> that contains the non zero entries of the upper triangle of 
		/// <paramref name="original"/>. In skyline format, some zero entries will be explicitly stored.
		/// </summary>
		/// <param name="original">The matrix that will be copied. It must be symmetric.</param>
		/// <param name="tolerance">
		/// The tolerance used to determine if an entry is zero. It will also be used to check, if <paramref name="original"/>
		/// is symmetric (<paramref name="original"/>[i, j] == <paramref name="original"/>[j, i]).
		/// </param>
		public static SkylineMatrix CreateFromArray(double[,] original, double tolerance = 1E-10)
		{
			if (!original.IsSymmetric(tolerance)) throw new ArgumentException("The original matrix must be symmetric.");

			// Indexing array
			var comparer = new ValueComparer(tolerance);
			int order = original.GetLength(0);
			var diagOffsets = new int[order + 1];
			//diagOffsets[0] = 0; // by default;
			for (int j = 0; j < order; ++j)
			{
				int colHeight = 0;
				for (int i = j - 1; i >= 0; --i)
				{
					if (!comparer.AreEqual(0.0, original[i, j])) colHeight = j - i;
				}
				diagOffsets[j + 1] = diagOffsets[j] + colHeight + 1;
			}

			// Values array
			int nnz = diagOffsets[order];
			var values = new double[nnz];
			for (int j = 0; j < order; ++j)
			{
				int colHeight = diagOffsets[j + 1] - diagOffsets[j] - 1;
				for (int t = 0; t <= colHeight; ++t)
				{
					int i = j - t; // row index of Aij
					int offsetAij = diagOffsets[j] + t;
					values[offsetAij] = original[i, j];
				}
			}

			return new SkylineMatrix(order, values, diagOffsets);
		}

		/// <summary>
		/// Initializes a new <see cref="SkylineMatrix"/> with the specified dimensions and the provided arrays 
		/// (<paramref name="values"/> and <paramref name="diagOffsets"/>) as its internal data.
		/// </summary>
		/// <param name="order">The number of rows/columns of the new matrix.</param>
		/// <param name="values">Contains the non zero superdiagonal entries of the matrix in column major order, starting from 
		///     the diagonal and going upwards.</param>
		/// <param name="diagOffsets">Contains the indices into <paramref name="values"/> of the diagonal entries of the matrix. 
		///     Its length is <paramref name="order"/> + 1, with the last entry being equal to nnz.</param>
		/// <param name="checkInput">If true, the provided arrays will be checked to make sure they are valid Skyline arrays, 
		///     which is safer. If false, no such check will take place, which is faster.</param>
		/// <param name="copyArrays">If true, the provided arrays will be copied and the new <see cref="SkylineMatrix"/> instance 
		///     will have references to the copies, which is safer. If false, the new matrix will have references to the 
		///     provided arrays themselves, which is faster.</param>
		public static SkylineMatrix CreateFromArrays(int order, double[] values, int[] diagOffsets,
			bool checkInput, bool copyArrays = false)
		{
			if (checkInput)
			{
				if (diagOffsets.Length != order + 1)
				{
					throw new ArgumentException("The length of the Skyline diagonal offsets array must be equal to the number of"
						+ " rows/columns + 1, but was " + diagOffsets.Length);
				}
				if (diagOffsets[diagOffsets.Length - 1] != values.Length)
				{
					throw new ArgumentException("The last entry of the Skyline diagonal offsets array must be equal to the number"
						+ " of non zero entries, but was " + diagOffsets[diagOffsets.Length - 1]);
				}
			}
			if (copyArrays)
			{
				double[] valuesCopy = new double[values.Length];
				Array.Copy(values, valuesCopy, values.Length);
				int[] diagOffsetsCopy = new int[diagOffsets.Length];
				Array.Copy(diagOffsets, diagOffsetsCopy, diagOffsets.Length);
				return new SkylineMatrix(order, valuesCopy, diagOffsetsCopy);
			}
			else return new SkylineMatrix(order, values, diagOffsets);
		}

		/// <summary>
		/// Initializes a new <see cref="SkylineMatrix"/> that contains the non zero entries of the upper triangle of 
		/// <paramref name="original"/>. In skyline format, some zero entries will be explicitly stored.
		/// </summary>
		/// <param name="original">The matrix that will be copied. It must be symmetric.</param>
		/// <param name="tolerance">
		/// The tolerance used to determine if an entry is zero. It will also be used to check, if <paramref name="original"/>
		/// is symmetric (<paramref name="original"/>[i, j] == <paramref name="original"/>[j, i]).
		/// </param>
		public static SkylineMatrix CreateFromMatrix(IIndexable2D original, double tolerance = 1E-10)
		{
			if (!original.IsSymmetric(tolerance)) throw new ArgumentException("The original matrix must be symmetric.");

			// Indexing array
			var comparer = new ValueComparer(tolerance);
			int order = original.NumColumns;
			var diagOffsets = new int[order + 1];
			//diagOffsets[0] = 0; // by default;
			for (int j = 0; j < order; ++j)
			{
				int colHeight = 0;
				for (int i = j - 1; i >= 0; --i)
				{
					if (!comparer.AreEqual(0.0, original[i, j])) colHeight = j - i;
				}
				diagOffsets[j + 1] = diagOffsets[j] + colHeight + 1;
			}

			// Values array
			int nnz = diagOffsets[order];
			var values = new double[nnz];
			for (int j = 0; j < order; ++j)
			{
				int colHeight = diagOffsets[j + 1] - diagOffsets[j] - 1;
				for (int t = 0; t <= colHeight; ++t)
				{
					int i = j - t; // row index of Aij
					int offsetAij = diagOffsets[j] + t;
					values[offsetAij] = original[i, j];
				}
			}

			return new SkylineMatrix(order, values, diagOffsets);
		}

		/// <summary>
		/// Initializes a new <see cref="SkylineMatrix"/> with the specified dimensions and the sparsity pattern defined by 
		/// <paramref name="diagOffsets"/>. The stored entries will initially be 0.
		/// </summary>
		/// <param name="order">The number of rows/columns of the new matrix.</param>
		/// <param name="diagOffsets">Contains the indices into <paramref name="values"/> of the diagonal entries of the matrix. 
		///     Its length is <paramref name="order"/> + 1, with the last entry being equal to nnz.</param>
		/// <param name="checkInput">If true, <paramref name="diagOffsets"/> will be checked to make sure it is a valid Skyline  
		///     array, which is safer. If false, no such check will take place, which is faster.</param>
		public static SkylineMatrix CreateZeroWithPattern(int order, int[] diagOffsets, bool checkInput)
		{
			if (checkInput)
			{
				if (diagOffsets.Length != order + 1)
				{
					throw new ArgumentException("The length of the Skyline diagonal offsets array must be equal to the number of"
						+ " rows/columns + 1, but was " + diagOffsets.Length);
				}
			}
			int nnz = diagOffsets[diagOffsets.Length] - 1;
			return new SkylineMatrix(order, new double[nnz], diagOffsets);
		}

		public override void AxpyIntoThis(SkylineMatrix otherMatrix, double otherCoefficient)
		{
			if (HasSameFormat(otherMatrix))
			{
				GlobalProvider.Blas.Daxpy(values.Length, otherCoefficient, otherMatrix.values, 0, 1, this.values, 0, 1);
			}
			else if (otherMatrix.values.Length == 0)
			{
				//TODO: I think this needs to throw an exception. When would this work? Both matrices must be zero.
				//		Otherwise the values array of the other matrix is overwritten, thus it is invalid. In any case, it is not
				//		the job of this matrix to operate on invalid matrices.
				Preconditions.CheckSameMatrixDimensions(this, otherMatrix);
				return; // The operation can be completed if the other matrix is empty.
			}
			else
			{
				Preconditions.CheckSameMatrixDimensions(this, otherMatrix);
				for (int j = 0; j < NumColumns; ++j)
				{
					// Check if the current column of this matrix is tall enough
					int thisDiagOffset = this.diagOffsets[j];
					int otherDiagOffset = otherMatrix.diagOffsets[j];
					int thisColTop = j - this.diagOffsets[j + 1] + thisDiagOffset + 1;
					int otherColTop = j - this.diagOffsets[j + 1] + otherDiagOffset + 1;
					if (thisColTop < otherColTop)
					{
						throw new SparsityPatternModifiedException(
							$"Column {j} of this matrix is shorter, which would result in overflow");
					}

					// Do the operation between the two columns. The column of this matrix is taller or equal to the other one.
					for (int i = j; i >= otherColTop; --i) // non zero entries of shortest=other column, including diagonal
					{
						this.values[thisDiagOffset + j - i] += otherCoefficient * otherMatrix.values[otherDiagOffset + j - i];
					}
					// Don't do anything to the non zero entries of the above the shortest=other column (this[i,j] += a*0)
				}
			}
		}

		public override SkylineMatrix CopyAsSameType(bool copyIndexingData)
		{
			var valuesCopy = new double[this.values.Length];
			Array.Copy(this.values, valuesCopy, this.values.Length);

			if (!copyIndexingData) return new SkylineMatrix(NumColumns, valuesCopy, this.diagOffsets);
			else
			{
				var diagOffsetsCopy = new int[this.diagOffsets.Length];
				Array.Copy(this.diagOffsets, diagOffsetsCopy, this.diagOffsets.Length);
				return new SkylineMatrix(NumColumns, valuesCopy, diagOffsetsCopy);
			}
		}

		/// <summary>
		/// Copies the entries of the matrix into a 2-dimensional array. The returned array has length(0) = <see cref="NumRows"/> 
		/// and length(1) = <see cref="NumColumns"/>.
		/// </summary>
		public double[,] CopyToArray2D()
		{
			double[,] array2D = new double[NumColumns, NumColumns];
			for (int j = 0; j < NumColumns; ++j)
			{
				int colOffset = diagOffsets[j];
				int columnTop = j - diagOffsets[j + 1] + colOffset + 1;
				array2D[j, j] = values[colOffset]; // diagonal entry
				for (int i = columnTop; i < j; ++i) // non zero entries stored above diagonal
				{
					double value = values[colOffset + j - i];
					array2D[j, i] = value;
					array2D[i, j] = value;
				}
			}
			return array2D;
		}

		public override Matrix CopyToFullMatrix()
		{
			Matrix fullMatrix = Matrix.CreateZero(this.NumColumns, this.NumColumns);
			for (int j = 0; j < NumColumns; ++j)
			{
				int colOffset = diagOffsets[j];
				int columnTop = j - diagOffsets[j + 1] + colOffset + 1;
				fullMatrix[j, j] = values[colOffset]; // diagonal entry
				for (int i = columnTop; i < j; ++i) // non zero entries stored above diagonal
				{
					double value = values[colOffset + j - i];
					fullMatrix[j, i] = value;
					fullMatrix[i, j] = value;
				}
			}
			return fullMatrix;
		}

		public int CountNonZeros() => values.Length;

		public override SkylineMatrix CreateZeroMatrixSame()
		{
			var resultValues = new double[values.Length];
			return new SkylineMatrix(NumColumns, resultValues, diagOffsets);
		}

		public override void DoEntrywiseIntoThis(SkylineMatrix otherMatrix, Func<double, double, double> binaryOperation)
		{
			if (HasSameFormat(otherMatrix))
			{
				for (int i = 0; i < values.Length; ++i)
				{
					this.values[i] = binaryOperation(this.values[i], otherMatrix.values[i]);
				}
			}
			else if (otherMatrix.values.Length == 0)
			{
				Preconditions.CheckSameMatrixDimensions(this, otherMatrix);
				return; // The operation can be completed if the other matrix is empty.
			}
			else
			{
				Preconditions.CheckSameMatrixDimensions(this, otherMatrix);
				for (int j = 0; j < NumColumns; ++j)
				{
					// Check if the current column of this matrix is tall enough
					int thisDiagOffset = this.diagOffsets[j];
					int otherDiagOffset = otherMatrix.diagOffsets[j];
					int thisColTop = j - this.diagOffsets[j + 1] + thisDiagOffset + 1;
					int otherColTop = j - this.diagOffsets[j + 1] + otherDiagOffset + 1;
					if (thisColTop < otherColTop)
					{
						throw new SparsityPatternModifiedException(
							$"Column {j} of this matrix is shorter, which would result in overflow");
					}

					// Do the operation between the two columns. The column of this matrix is taller or equal to the other one.
					for (int i = j; i >= otherColTop; --i) // non zero entries of shortest column, including diagonal
					{
						int thisIndex = thisDiagOffset + j - i;
						this.values[thisIndex] = binaryOperation(
							this.values[thisIndex], otherMatrix.values[otherDiagOffset + j - i]);
					}

					for (int i = otherColTop - 1; i >= thisColTop; --i) // non zero entries of the above the shortest column
					{
						int thisIndex = thisDiagOffset + j - i;
						this.values[thisIndex] = binaryOperation(
							this.values[thisIndex], otherMatrix.values[otherDiagOffset + j - i]);
					}
				}
			}
		}

		public IEnumerable<(int row, int col, double value)> EnumerateNonZeros()
		{
			for (int j = 0; j < NumColumns; ++j)
			{
				int colOffset = diagOffsets[j];
				int columnTop = j - diagOffsets[j + 1] + colOffset + 1;
				yield return (j, j, values[colOffset]); // diagonal entry
				for (int i = columnTop; i < j; ++i)
				{
					double value = values[colOffset + j - i];
					yield return (i, j, value);
					yield return (j, i, value);
				}
			}
		}

		public override bool Equals(IIndexable2D other, double tolerance = 1E-13)
		{
			if ((this.NumRows != other.NumRows) || (this.NumColumns != other.NumColumns)) return false;
			var comparer = new ValueComparer(1e-13);
			for (int j = 0; j < NumColumns; ++j)
			{
				int colOffset = diagOffsets[j];
				int columnTop = j - diagOffsets[j + 1] + colOffset + 1;
				for (int i = 0; i < columnTop; ++i) // zero entries above stored column
				{
					if (!(comparer.AreEqual(0.0, other[i, j]) && comparer.AreEqual(0.0, other[j, i]))) return false;
				}
				for (int i = columnTop; i < j; ++i) // non zero entries of column, excluding diafonal
				{
					double value = values[colOffset + j - i];
					if (!(comparer.AreEqual(value, other[i, j]) && comparer.AreEqual(value, other[j, i]))) return false;
				}
				if (!comparer.AreEqual(values[colOffset], other[j, j])) return false; // non zero diagonal entry
			}
			return true; // At this point all entries have been checked and are equal
		}

		/// <summary>
		/// Calculate the Cholesky factorization. The matrix must be positive definite, otherwise an
		/// <see cref="IndefiniteMatrixException"/> will be thrown. If <paramref name="inPlace"/> is set to true, this object 
		/// must not be used again, otherwise a <see cref="NullReferenceException"/> will be thrown.
		/// </summary>
		/// <param name="inPlace">
		/// False, to copy the internal non zero entries before factorization. True, to overwrite them with the factorized data, 
		/// thus saving memory and time. However, that will make this object unusable, so you MUST NOT call any other members 
		/// afterwards.
		/// </param>
		/// <param name="pivotTolerance">
		/// If a diagonal entry is closer to zero than this tolerance, an <see cref="IndefiniteMatrixException"/> exception will
		/// be thrown.
		/// </param>
		/// <exception cref="IndefiniteMatrixException">Thrown if the matrix is not positive definite.</exception>
		/// <exception cref="NullReferenceException">
		/// Thrown if a member of his instance is accessed after this method is called.
		/// </exception>"
		public CholeskySkyline FactorCholesky(bool inPlace, double tolerance = CholeskySkyline.PivotTolerance)
		{
			if (inPlace)
			{
				var factor = CholeskySkyline.Factorize(NumColumns, values, diagOffsets, tolerance);
				// Set the skyline arrays to null to force NullReferenceException if they are accessed again.
				// TODO: perhaps there is a better way to handle this.
				values = null;
				diagOffsets = null;
				isOverwritten = true;
				return factor;
			}
			else
			{
				double[] valuesCopy = new double[values.Length];
				Array.Copy(values, valuesCopy, values.Length);
				return CholeskySkyline.Factorize(NumColumns, valuesCopy, diagOffsets, tolerance);
			}
		}

		/// <summary>
		/// Calculate the LDL factorization. The matrix must be invertible, otherwise a <see cref="SingularMatrixException"/> 
		/// will be thrown. This method succeeds for all symmetric positive definite matrices, but not for all symmetric ones. 
		/// If <paramref name="inPlace"/> is set to true, this object must not be used again, otherwise a 
		/// <see cref="NullReferenceException"/> will be thrown.
		/// </summary>
		/// <param name="inPlace">
		/// False, to copy the internal non zero entries before factorization. True, to overwrite them with the factorized data, 
		/// thus saving memory and time. However, that will make this object unusable, so you MUST NOT call any other members 
		/// afterwards.
		/// </param>
		/// <param name="tolerance">
		/// If a diagonal entry is closer to zero than this tolerance, an <see cref="SingularMatrixException"/> exception will 
		/// be thrown.
		/// </param>
		/// <exception cref="SingularMatrixException">Thrown if the matrix is not invertible.</exception>
		/// <exception cref="NullReferenceException">
		/// Thrown if a member of his instance is accessed after this method is called.
		/// </exception>"
		public LdlSkyline FactorLdl(bool inPlace, double tolerance = LdlSkyline.PivotTolerance)
		{
			if (inPlace)
			{
				var factor = LdlSkyline.Factorize(NumColumns, values, diagOffsets, tolerance);
				// Set the skyline arrays to null to force NullReferenceException if they are accessed again.
				// TODO: perhaps there is a better way to handle this.
				values = null;
				diagOffsets = null;
				isOverwritten = true;
				return factor;
			}
			else
			{
				double[] valuesCopy = new double[values.Length];
				Array.Copy(values, valuesCopy, values.Length);
				return LdlSkyline.Factorize(NumColumns, valuesCopy, diagOffsets, tolerance);
			}
		}

		/// <summary>
		/// Applies the Cholesky factorization to the independent columns of a symmetric positive semi-definite matrix,
		/// sets the dependent ones equal to columns of the identity matrix and return the nullspace of the matrix. Requires 
		/// extra memory for the basis vectors of the nullspace. If <paramref name="inPlace"/> is set to true, this object 
		/// must not be used again, otherwise a <see cref="NullReferenceException"/> will be thrown.
		/// </summary>
		/// <param name="inPlace">
		/// False, to copy the internal non zero entries before factorization. True, to overwrite them with the factorized data, 
		/// thus saving memory and time. However, that will make this object unusable, so you MUST NOT call any other members 
		/// afterwards.
		/// </param>
		/// <param name="pivotTolerance">
		/// If a diagonal entry is &lt;= <paramref name="pivotTolerance"/> it means that the corresponding column is dependent 
		/// on the rest. The Cholesky factorization only applies to independent column, while dependent ones are used to compute
		/// the nullspace. Therefore it is important to select a tolerance that will identify small pivots that result from 
		/// singularity, but not from ill-conditioning.
		/// </param>
		/// <exception cref="IndefiniteMatrixException">Thrown if the matrix is not positive definite.</exception>
		/// <exception cref="NullReferenceException">
		/// Thrown if a member of his instance is accessed after this method is called.
		/// </exception>"
		public SemidefiniteCholeskySkyline FactorSemidefiniteCholesky(bool inPlace,
			double pivotTolerance = SemidefiniteCholeskySkyline.PivotTolerance)
		{
			if (inPlace)
			{
				var factor = SemidefiniteCholeskySkyline.Factorize(NumColumns, values, diagOffsets, pivotTolerance);
				// Set the skyline arrays to null to force NullReferenceException if they are accessed again.
				// TODO: perhaps there is a better way to handle this.
				values = null;
				diagOffsets = null;
				isOverwritten = true;
				return factor;
			}
			else
			{
				double[] valuesCopy = new double[values.Length];
				Array.Copy(values, valuesCopy, values.Length);
				return SemidefiniteCholeskySkyline.Factorize(NumColumns, valuesCopy, diagOffsets, pivotTolerance);
			}
		}

		/// <summary>
		/// Applies the LDL factorization to the independent columns of a symmetric positive semi-definite matrix,
		/// sets the dependent ones equal to columns of the identity matrix and return the nullspace of the matrix. Requires 
		/// extra memory for the basis vectors of the nullspace. If <paramref name="inPlace"/> is set to true, this object 
		/// must not be used again, otherwise a <see cref="NullReferenceException"/> will be thrown.
		/// </summary>
		/// <param name="inPlace">
		/// False, to copy the internal non zero entries before factorization. True, to overwrite them with the factorized data, 
		/// thus saving memory and time. However, that will make this object unusable, so you MUST NOT call any other members 
		/// afterwards.
		/// </param>
		/// <param name="pivotTolerance">
		/// If a diagonal entry is &lt;= <paramref name="pivotTolerance"/> it means that the corresponding column is dependent 
		/// on the rest. The Cholesky factorization only applies to independent column, while dependent ones are used to compute
		/// the nullspace. Therefore it is important to select a tolerance that will identify small pivots that result from 
		/// singularity, but not from ill-conditioning.
		/// </param>
		/// <exception cref="IndefiniteMatrixException">Thrown if the matrix is not positive definite.</exception>
		/// <exception cref="NullReferenceException">
		/// Thrown if a member of his instance is accessed after this method is called.
		/// </exception>"
		public SemidefiniteLdlSkyline FactorSemidefiniteLdl(bool inPlace,
			double pivotTolerance = SemidefiniteLdlSkyline.PivotTolerance)
		{
			if (inPlace)
			{
				var factor = SemidefiniteLdlSkyline.Factorize(NumColumns, values, diagOffsets, pivotTolerance);
				// Set the skyline arrays to null to force NullReferenceException if they are accessed again.
				// TODO: perhaps there is a better way to handle this.
				values = null;
				diagOffsets = null;
				isOverwritten = true;
				return factor;
			}
			else
			{
				double[] valuesCopy = new double[values.Length];
				Array.Copy(values, valuesCopy, values.Length);
				return SemidefiniteLdlSkyline.Factorize(NumColumns, valuesCopy, diagOffsets, pivotTolerance);
			}
		}

		public override Vector GetColumn(int colIndex)
		{
			if (isOverwritten) throw new MatrixDataOverwrittenException();
			Preconditions.CheckIndexCol(this, colIndex);
			return Vector.CreateFromArray(SkylineSlicing.GetColumn(values, diagOffsets, colIndex));
		}

		public override double[] GetDiagonalAsArray()
		{
			if (isOverwritten) throw new MatrixDataOverwrittenException();
			return SkylineSlicing.GetDiagonal(values, diagOffsets);
		}

		public override Vector GetRow(int rowIndex) => GetColumn(rowIndex);

		public SparseFormat GetSparseFormat()
		{
			if (isOverwritten) throw new MatrixDataOverwrittenException();
			var format = new SparseFormat();
			format.RawValuesTitle = "Values";
			format.RawValuesArray = values;
			format.RawIndexArrays.Add("Diagonal offsets", diagOffsets);
			return format;
		}

		public override IMatrix GetSubmatrix(
			int rowStartInclusive, int rowEndExclusive, int colStartInclusive, int colEndExclusive)
		{
			if (isOverwritten) throw new MatrixDataOverwrittenException();
			int[] rowIndices = Enumerable.Range(rowStartInclusive, rowEndExclusive - rowStartInclusive).ToArray();
			int[] colIndices = Enumerable.Range(colStartInclusive, colEndExclusive - colStartInclusive).ToArray();
			return GetSubmatrix(rowIndices, colIndices);
		}

		public CscMatrix GetSubmatrixCsc(int[] rowIndices, int[] colIndices)
		{
			if (isOverwritten) throw new MatrixDataOverwrittenException();
			return SkylineSlicing.GetSubmatrixCsc(values, diagOffsets, rowIndices, colIndices);
		}

		public Matrix GetSubmatrixSymmetricFull(int[] indices)
		{
			if (isOverwritten) throw new MatrixDataOverwrittenException();
			return SkylineSlicing.GetSubmatrixSymmetricFull(values, diagOffsets, indices);
			//// I am not sure that the above is faster than: 
			//return DenseStrategies.GetSubmatrix(this, indices, indices);
		}

		public SymmetricMatrix GetSubmatrixSymmetricPacked(int[] indices)
		{
			//TODO: perhaps this can be combined with the CSC and full version to get all 2 submatrices needed for 
			//      Schur complements more efficiently.
			if (isOverwritten) throw new MatrixDataOverwrittenException();
			return SkylineSlicing.GetSubmatrixSymmetricPacked(values, diagOffsets, indices);
		}

		public SparsityPatternSymmetric GetSubmatrixSymmetricPattern(int[] indices)
		{
			//TODO: perhaps this can be combined with the CSC and full version to get all 2 submatrices needed for 
			//      Schur complements more efficiently.
			if (isOverwritten) throw new MatrixDataOverwrittenException();
			return SkylineSlicing.GetSubmatrixSymmetricPattern(values, diagOffsets, indices);
		}

		public SkylineMatrix GetSubmatrixSymmetricSkyline(int[] indices)
		{
			//TODO: perhaps this can be combined with the CSC and full version to get all 2 submatrices needed for 
			//      Schur complements more efficiently.
			if (isOverwritten) throw new MatrixDataOverwrittenException();
			return SkylineSlicing.GetSubmatrixSymmetricSkyline(values, diagOffsets, indices);
		}

		public override bool HasSameFormat(SkylineMatrix otherMatrix) => this.diagOffsets == otherMatrix.diagOffsets;

		public override void LinearCombinationIntoThis(double thisCoefficient, SkylineMatrix otherMatrix, double otherCoefficient)
		{
			if (HasSameFormat(otherMatrix))
			{
				GlobalProvider.Blas.Daxpby(
					RawValues.Length, otherCoefficient, otherMatrix.RawValues, 0, 1, thisCoefficient, this.RawValues, 0, 1);
			}
			else if (otherMatrix.RawValues.Length == 0)
			{
				Preconditions.CheckSameMatrixDimensions(this, otherMatrix);
				return; // The operation can be completed if the other matrix is empty.
			}
			else
			{
				Preconditions.CheckSameMatrixDimensions(this, otherMatrix);
				for (int j = 0; j < NumColumns; ++j)
				{
					// Check if the current column of this matrix is tall enough
					int thisDiagOffset = this.diagOffsets[j];
					int otherDiagOffset = otherMatrix.diagOffsets[j];
					int thisColTop = j - this.diagOffsets[j + 1] + thisDiagOffset + 1;
					int otherColTop = j - this.diagOffsets[j + 1] + otherDiagOffset + 1;
					if (thisColTop < otherColTop)
					{
						throw new SparsityPatternModifiedException(
							$"Column {j} of this matrix is shorter, which would result in overflow");
					}

					// Do the operation between the two columns. The column of this matrix is taller or equal to the other one.
					for (int i = j; i >= otherColTop; --i) // non zero entries of shortest column, including diagonal
					{
						int thisIndex = thisDiagOffset + j - i;
						this.values[thisIndex] = thisCoefficient * this.values[thisIndex]
							+ otherCoefficient * otherMatrix.values[otherDiagOffset + j - i];
					}

					for (int i = otherColTop - 1; i >= thisColTop; --i) // non zero entries of the above the shortest column
					{
						int thisIndex = thisDiagOffset + j - i;
						this.values[thisIndex] = thisCoefficient * this.values[thisIndex];
					}
				}
			}
		}

		/// <summary>
		/// See <see cref="IReadOnlyMatrix.Multiply(IReadOnlyVector, bool)"/>.
		/// </summary>
		/// <remarks>
		/// <paramref name="transposeThis"/> does not affect the result, as a <see cref="SkylineMatrix"/> is symmetric.
		/// </remarks>
		public override IVector Multiply(IReadOnlyVector vector, bool transposeThis = false)
		{
			if (vector is Vector casted) return Multiply(casted);
			else throw new NotImplementedException();
		}

		/// <summary>
		/// Performs the matrix-vector multiplication: this * <paramref name="vector"/> = <paramref name="vector"/> * this.
		/// </summary>
		/// <param name="vector">A vector with <see cref="IIndexable1D.Length"/> being equal to 
		///     this.<see cref="NumColumns"/>.</param>
		/// <exception cref="NonMatchingDimensionsException">Thrown if the <see cref="IIndexable1D.Length"/> of
		///     <paramref name="vector"/> is different than the <see cref="NumColumns"/> of this.</exception>
		public Vector Multiply(Vector vector)
		{
			//TODO: this performs redundant dimension checks
			var result = Vector.CreateZero(NumColumns);
			MultiplyIntoResult(vector, result);
			return result;
		}

		/// <summary>
		/// See <see cref="IReadOnlyMatrix.MultiplyIntoResult(IReadOnlyVector, IVector, bool)"/>.
		/// </summary>
		public override void MultiplyIntoResult(IReadOnlyVector lhsVector, IVector rhsVector, bool transposeThis = false)
		{
			if (this.values.Length == 0)
			{
				return;
			}

			if ((lhsVector is Vector lhsDense) && (rhsVector is Vector rhsDense))
			{
				MultiplyIntoResult(lhsDense, rhsDense);
			}
			else throw new NotImplementedException();
		}

		/// <summary>
		/// Performs the matrix-vector multiplication: <paramref name="rhsVector"/> = this * <paramref name="vector"/>.
		/// To multiply this * columnVector, set <paramref name="transposeThis"/> to false.
		/// To multiply rowVector * this, set <paramref name="transposeThis"/> to true.
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
			Preconditions.CheckMultiplicationDimensions(NumColumns, lhsVector.Length);
			Preconditions.CheckSystemSolutionDimensions(NumRows, rhsVector.Length);
			if (this.values.Length == 0)
			{
				return;
			}

			ManagedSparseBlasProvider.UniqueInstance.Dskymv(
				NumColumns, values, diagOffsets, lhsVector.RawData, rhsVector.RawData);
		}

		public override double Reduce(double identityValue, ProcessEntry processEntry, ProcessZeros processZeros, Finalize finalize)
		{
			double aggregator = identityValue;
			int n = NumColumns;
			int nnz = values.Length;
			for (int j = 0; j < n; ++j)
			{
				int colStart = diagOffsets[j];
				int colEnd = diagOffsets[j + 1];
				processEntry(values[colStart], aggregator); // process diagonal entry once
				for (int t = colStart + 1; t < colEnd; ++t) // process off-diagonal entries twice (once for each triangle)
				{
					processEntry(values[t], aggregator);
					processEntry(values[t], aggregator);
				}
			}

			aggregator = processZeros(n * n - nnz, aggregator);
			return finalize(aggregator);
		}

		public override IMatrix Transpose() => CopyAsSameType(true);

		/// <summary>
		/// Perhaps this should be manually inlined. Testing needed.
		/// </summary>
		/// <param name="rowIdx"></param>
		/// <param name="colIdx"></param>
		[MethodImpl(MethodImplOptions.AggressiveInlining)]
		internal static void ProcessIndices(ref int rowIdx, ref int colIdx)
		{
			if (rowIdx > colIdx)
			{
				int swap = rowIdx;
				rowIdx = colIdx;
				colIdx = swap;
			}
		}
	}
}
