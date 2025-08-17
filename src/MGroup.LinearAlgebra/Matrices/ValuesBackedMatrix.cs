namespace MGroup.LinearAlgebra.Matrices
{
	using System;
	using System.Collections.Generic;
	using System.Text;

	using MGroup.LinearAlgebra.Commons;
	using MGroup.LinearAlgebra.Exceptions;
	using MGroup.LinearAlgebra.Reduction;

	using static MGroup.LinearAlgebra.LibrarySettings;

	[Serializable]
	public abstract class ValuesBackedMatrix<TMatrix> : DefaultMatrix
		where TMatrix : ValuesBackedMatrix<TMatrix>
	{
		/// <summary>
		/// The internal array that stores the non-zero entries of the matrix. The non-zero entries of each row are consecutive.
		/// Its length is equal to the number of non-zero entries.
		/// It should only be used for passing the raw array to linear algebra libraries.
		/// </summary>
		public abstract double[] RawValues { get; }

		/// <summary>
		/// Copies the entries of this matrix.
		/// </summary>
		/// <param name="copyIndexingData">
		/// If true, all data of this object will be copied. If false, only the array containing the values of the stored
		/// matrix entries will be copied. The new matrix will reference the same indexing arrays as this one.
		/// </param>
		/// <returns>A matrix with the same type as this.</returns>
		public abstract TMatrix CopyAsSameType(bool copyIndexingData);

		/// <summary>
		/// Initializes a new instance of the same type as this matrix, with the exact same storage format and zero entries.
		/// </summary>
		/// <returns>A matrix with the same sparsity pattern as this instance, but all its explitly stored values are 0.</returns>
		public abstract TMatrix CreateZeroMatrixSame();

		/// <summary>
		/// Returns true if the only difference betweens this matrix and <paramref name="other"/> is their values.
		/// </summary>
		/// <param name="other">The matrix to compare.</param>
		/// <returns>True if the matrices have the same format. False otherwise.</returns>
		public abstract bool HasSameFormat(TMatrix other);

		public override IMatrix Axpy(IMatrixView otherMatrix, double otherCoefficient)
		{
			if (otherMatrix is TMatrix casted && HasSameFormat(casted))
			{
				TMatrix result = CopyAsSameType(false);
				result.AxpyIntoThis(casted, otherCoefficient);
				return result;
			}
			else
			{
				// All entries must be processed. TODO: optimizations may be possible (e.g. only access the nnz in this matrix)
				return base.Axpy(this, otherCoefficient);
			}
		}

		public override void AxpyIntoThis(IMatrixView otherMatrix, double otherCoefficient)
		{
			if (otherMatrix is TMatrix casted)
			{
				AxpyIntoThis(casted, otherCoefficient);
			}
			else
			{
				throw new SparsityPatternModifiedException(
					"This operation is allowed only if the other matrix has the same sparsity pattern");
			}
		}

		/// <summary>
		/// For each pair of entries (i, j) performs the operation:
		/// this[i, j] = <paramref name="otherCoefficient"/> * <paramref name="otherMatrix"/>[i, j] + this[i, j].
		/// The resulting matrix overwrites the entries of this <typeparamref name="TMatrix"/> instance.
		/// </summary>
		/// <param name="otherMatrix">
		/// A matrix with the same indexing arrays as this <typeparamref name="TMatrix"/> instance.
		/// </param>
		/// <param name="otherCoefficient">
		/// A scalar that multiplies each entry of <paramref name="otherMatrix"/>.
		/// </param>
		/// <exception cref="SparsityPatternModifiedException">
		/// Thrown if <paramref name="otherMatrix"/> has different indexing arrays than this instance.
		/// </exception>
		public virtual void AxpyIntoThis(TMatrix otherMatrix, double otherCoefficient)
		{
			if (HasSameFormat(otherMatrix))
			{
				GlobalProvider.Blas.Daxpy(RawValues.Length, otherCoefficient, otherMatrix.RawValues, 0, 1, this.RawValues, 0, 1);
			}
			else if (otherMatrix.RawValues.Length == 0)
			{
				//TODO: I think this needs to throw an exception. When would this work? Both matrices must be zero.
				//		Otherwise the values array of the other matrix is overwritten, thus it is invalid. In any case, it is not
				//		the job of this matrix to operate on invalid matrices.
				Preconditions.CheckSameMatrixDimensions(this, otherMatrix);
				return; // The operation can be completed if the other matrix is empty.
			}
			else
			{
				throw new SparsityPatternModifiedException(
					"This operation is allowed only if the other matrix has the same sparsity pattern");
			}
		}

		/// <summary>
		/// For each pair of entries (i, j) performs the operation:
		/// result[i, j] = <paramref name="otherCoefficient"/> * <paramref name="otherMatrix"/>[i, j] + this[i, j].
		/// The resulting matrix is written to a new <typeparamref name="TMatrix"/> and then returned.
		/// </summary>
		/// <param name="otherMatrix">
		/// A matrix with the same number of rows and columns and the same type as this matrix instance.
		/// </param>
		/// <param name="otherCoefficient">A scalar that multiplies each entry of <paramref name="otherMatrix"/>.</param>
		/// <exception cref="SparsityPatternModifiedException">
		/// Thrown if <paramref name="otherMatrix"/> has different indexing arrays than this instance.
		/// </exception>
		/// <returns>A matrix with the same type as this.</returns>
		public virtual TMatrix AxpySameFormat(TMatrix otherMatrix, double otherCoefficient)
		{
			if (HasSameFormat(otherMatrix))
			{
				TMatrix result = CopyAsSameType(false);
				result.AxpyIntoThis(otherMatrix, otherCoefficient);
				return result;
			}
			else
			{
				// Conceptually it is not wrong to do this, even if the indexers are different, but how would I implement it?
				throw new SparsityPatternModifiedException(
					"This operation is allowed only if the other matrix has the same sparsity pattern");
			}
		}

		public override void Clear() => Array.Clear(RawValues, 0, RawValues.Length);

		public override IMatrix Copy(bool copyIndexingData = false) => CopyAsSameType(copyIndexingData);

		public override IMatrix CreateZeroMatrixWithSameFormat() => CreateZeroMatrixSame();

		public override IMatrix DoEntrywise(IMatrixView other, Func<double, double, double> binaryOperation)
		{
			if (other is TMatrix casted && HasSameFormat(casted))
			{
				TMatrix result = CopyAsSameType(false);
				result.DoEntrywiseIntoThis(casted, binaryOperation);
				return result;
			}
			else
			{
				// All entries must be processed. TODO: optimizations may be possible (e.g. only access the nnz in this matrix)
				return base.DoEntrywise(other, binaryOperation);
			}
		}

		public override void DoEntrywiseIntoThis(IMatrixView other, Func<double, double, double> binaryOperation)
		{
			if (other is TMatrix casted)
			{
				DoEntrywiseIntoThis(casted, binaryOperation);
			}
			else
			{
				throw new SparsityPatternModifiedException(
					"This operation is allowed only if the other matrix has the same sparsity pattern");
			}
		}

		/// <summary>
		/// For each pair of entries (i, j) performs the operation:
		/// this[i, j] = <paramref name="binaryOperation"/>(this[i,j], <paramref name="otherMatrix"/>[i,j]).
		/// The resulting matrix overwrites the entries of this.
		/// </summary>
		/// <param name="otherMatrix">
		/// A matrix with the same indexing arrays as this <typeparamref name="TMatrix"/> instance.
		/// </param>
		/// <param name="binaryOperation">A method that takes 2 arguments and returns 1 result.</param>
		/// <exception cref="SparsityPatternModifiedException">
		/// Thrown if <paramref name="otherMatrix"/> has different indexing arrays than this instance.
		/// </exception>
		public virtual void DoEntrywiseIntoThis(TMatrix otherMatrix, Func<double, double, double> binaryOperation)
		{
			if (HasSameFormat(otherMatrix))
			{
				for (int i = 0; i < RawValues.Length; ++i)
				{
					this.RawValues[i] = binaryOperation(this.RawValues[i], otherMatrix.RawValues[i]);
				}
			}
			else if (otherMatrix.RawValues.Length == 0)
			{
				Preconditions.CheckSameMatrixDimensions(this, otherMatrix);
				return; // The operation can be completed if the other matrix is empty.
			}
			else
			{
				throw new SparsityPatternModifiedException(
					"This operation is allowed only if the other matrix has the same sparsity pattern");
			}
		}

		/// <summary>
		/// For each pair of entries (i, j) performs the operation:
		/// this[i, j] = <paramref name="binaryOperation"/>(this[i,j], <paramref name="otherMatrix"/>[i,j]).
		/// The resulting matrix is written to a new <typeparamref name="TMatrix"/> and then returned.
		/// </summary>
		/// <param name="otherMatrix">
		/// A matrix with the same number of rows and columns and the same type as this matrix instance.
		/// </param>
		/// <param name="binaryOperation">A method that takes 2 arguments and returns 1 result.</param>
		/// <exception cref="SparsityPatternModifiedException">
		/// Thrown if <paramref name="otherMatrix"/> has different indexing arrays than this instance.
		/// </exception>
		/// <returns>A matrix with the same type as this.</returns>
		public virtual TMatrix DoEntrywiseSameFormat(TMatrix otherMatrix, Func<double, double, double> binaryOperation)
		{
			if (HasSameFormat(otherMatrix))
			{
				TMatrix result = CopyAsSameType(false);
				result.DoEntrywiseIntoThis(otherMatrix, binaryOperation);
				return result;
			}
			else
			{
				// Conceptually it is not wrong to do this, even if the indexers are different, but how would I implement it?
				throw new SparsityPatternModifiedException(
					"This operation is allowed only if the other matrix has the same sparsity pattern");
			}
		}

		public override IMatrix DoToAllEntries(Func<double, double> unaryOperation)
		{
			if (unaryOperation(0.0) == 0.0)
			{
				TMatrix result = CopyAsSameType(false);
				result.DoToAllEntriesIntoThis(unaryOperation);
				return result;
			}
			else
			{
				Matrix result = CopyToFullMatrix();
				result.DoToAllEntriesIntoThis(unaryOperation);
				return result;
			}
		}

		public override void DoToAllEntriesIntoThis(Func<double, double> unaryOperation)
		{
			if (unaryOperation(0.0) == 0.0)
			{
				for (int i = 0; i < RawValues.Length; ++i)
				{
					RawValues[i] = unaryOperation(RawValues[i]);
				}
			}
			else
			{
				throw new SparsityPatternModifiedException("This operation would change the sparsity pattern");
			}
		}

		/// <summary>
		/// For each pair of nonzero entries (i, j) performs the operation:
		/// result[i, j] = <paramref name="unaryOperation"/>(this[i, j]).
		/// The resulting matrix is written to a new <typeparamref name="TMatrix"/> and then returned.
		/// </summary>
		/// <param name="unaryOperation">A method that takes 1 argument and returns 1 result.</param>
		/// <returns>A matrix with the same type as this.</returns>
		public TMatrix DoToAllEntriesSameFormat(Func<double, double> unaryOperation)
		{
			if (unaryOperation(0.0) == 0.0)
			{
				TMatrix result = CopyAsSameType(false);
				result.DoToAllEntriesIntoThis(unaryOperation);
				return result;
			}
			else
			{
				throw new SparsityPatternModifiedException("This operation would change the sparsity pattern");
			}
		}

		public override bool HasSameFormat(IMatrixView otherMatrix)
		{
			if (otherMatrix is TMatrix casted)
			{
				return HasSameFormat(casted);
			}

			return false;
		}

		public override IMatrix LinearCombination(double thisCoefficient, IMatrixView otherMatrix, double otherCoefficient)
		{
			if (otherMatrix is TMatrix casted && HasSameFormat(casted))
			{
				TMatrix result = CopyAsSameType(false);
				result.LinearCombinationIntoThis(thisCoefficient, casted, otherCoefficient);
				return result;
			}
			else
			{
				// All entries must be processed. TODO: optimizations may be possible (e.g. only access the nnz in this matrix)
				return base.Axpy(this, otherCoefficient);
			}
		}

		public override void LinearCombinationIntoThis(double thisCoefficient, IMatrixView otherMatrix, double otherCoefficient)
		{
			if (otherMatrix is TMatrix casted)
			{
				LinearCombinationIntoThis(thisCoefficient, casted, otherCoefficient);
			}
			else
			{
				throw new SparsityPatternModifiedException(
					"This operation is allowed only if the other matrix has the same sparsity pattern");
			}
		}

		/// <summary>
		/// For each pair of entries (i, j) performs the operation: this[i, j] = <paramref name="thisCoefficient"/> * this[i, j]
		///  + <paramref name="otherCoefficient"/> * <paramref name="otherMatrix"/>[i, j].
		/// The resulting matrix overwrites the entries of this <typeparamref name="TMatrix"/> instance.
		/// </summary>
		/// <param name="thisCoefficient">A scalar that multiplies each entry of this matrix.</param>
		/// <param name="otherMatrix">
		/// A matrix with the same indexing arrays as this <typeparamref name="TMatrix"/> instance.
		/// </param>
		/// <param name="otherCoefficient">
		/// A scalar that multiplies each entry of <paramref name="otherMatrix"/>.
		/// </param>
		/// <exception cref="SparsityPatternModifiedException">
		/// Thrown if <paramref name="otherMatrix"/> has different indexing arrays than this instance.
		/// </exception>
		public virtual void LinearCombinationIntoThis(double thisCoefficient, TMatrix otherMatrix, double otherCoefficient)
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
				throw new SparsityPatternModifiedException(
					"This operation is allowed only if the other matrix has the same sparsity pattern");
			}
		}

		/// <summary>
		/// For each pair of entries (i, j) performs the operation: this[i, j] = <paramref name="thisCoefficient"/> * this[i, j]
		///  + <paramref name="otherCoefficient"/> * <paramref name="otherMatrix"/>[i, j].
		/// The resulting matrix is written to a new <typeparamref name="TMatrix"/> and then returned.
		/// </summary>
		/// <param name="thisCoefficient">A scalar that multiplies each entry of this matrix.</param>
		/// <param name="otherMatrix">
		/// A matrix with the same number of rows and columns and the same type as this matrix instance.
		/// </param>
		/// <param name="otherCoefficient">A scalar that multiplies each entry of <paramref name="otherMatrix"/>.</param>
		/// <exception cref="SparsityPatternModifiedException">
		/// Thrown if <paramref name="otherMatrix"/> has different indexing arrays than this instance.
		/// </exception>
		/// <returns>A matrix with the same type as this.</returns>
		public virtual TMatrix LinearCombinationSameFormat(double thisCoefficient, TMatrix otherMatrix, double otherCoefficient)
		{
			if (HasSameFormat(otherMatrix))
			{
				TMatrix result = CopyAsSameType(false);
				result.LinearCombinationIntoThis(thisCoefficient, otherMatrix, otherCoefficient);
				return result;
			}
			else
			{
				// Conceptually it is not wrong to do this, even if the indexers are different, but how would I implement it?
				throw new SparsityPatternModifiedException(
					"This operation is allowed only if the other matrix has the same sparsity pattern");
			}
		}

		public override IMatrix Scale(double scalar) => ScaleSameFormat(scalar);

		public override void ScaleIntoThis(double scalar) => GlobalProvider.Blas.Dscal(RawValues.Length, scalar, RawValues, 0, 1);

		/// <summary>
		/// For each pair of nonzero entries (i, j) performs the operation: result[i, j] = <paramref name="scalar"/> * this[i, j].
		/// The resulting matrix is written to a new <typeparamref name="TMatrix"/> and then returned.
		/// </summary>
		/// <param name="scalar">A scalar that multiplies each entry of this matrix.</param>
		public TMatrix ScaleSameFormat(double scalar)
		{
			TMatrix result = CopyAsSameType(false);
			result.ScaleIntoThis(scalar);
			return result;
		}

		protected double ReduceNonSymmetric(
			double identityValue, ProcessEntry processEntry, ProcessZeros processZeros, Finalize finalize)
		{
			double aggregator = identityValue;
			int nnz = RawValues.Length;
			for (int i = 0; i < nnz; ++i)
			{
				aggregator = processEntry(RawValues[i], aggregator);
			}

			aggregator = processZeros(NumRows * NumColumns - nnz, aggregator);
			return finalize(aggregator);
		}
	}
}
