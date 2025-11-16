using System;
using System.Collections.Generic;
using System.Linq;

using CSparse.Double.Factorization;

using MGroup.LinearAlgebra.Commons;
using MGroup.LinearAlgebra.Exceptions;
using MGroup.LinearAlgebra.Reduction;

using static MGroup.LinearAlgebra.LibrarySettings;

namespace MGroup.LinearAlgebra.Vectors
{
	/// <summary>
	/// A vector that only stores non-zero entries. Some zero entries can also be stored but they are  non-structural zeros 
	/// and thus handled as non-zero entries.
	/// </summary>
	public sealed class SparseVector : DefaultVector
	{
		private readonly double[] values;

		/// <summary>
		/// They must be sorted.
		/// </summary>
		private readonly int[] indices;

		private SparseVector(int length, double[] values, int[] indices)
		{
			this.Length = length;
			this.values = values;
			this.indices = indices;
		}

		public override int Length { get; }

		/// <summary>
		/// The internal array that stores the indices of the non-zero entries of the vector. 
		/// Its length is equal to the number of non-zero entries.
		/// It should only be used for passing the raw array to linear algebra libraries.
		/// </summary>
		public int[] RawIndices => indices;

		/// <summary>
		/// The internal array that stores the values of the non-zero entries of the vector.
		/// Its length is equal to the number of non-zero entries.
		/// It should only be used for passing the raw array to linear algebra libraries.
		/// </summary>
		public double[] RawValues => values;

		public override double this[int index]
		{
			get
			{
				int sparseIdx = FindSparseIndexOf(index);
				return (sparseIdx < 0) ? 0.0 : values[sparseIdx];
			}

			set
			{
				int sparseIdx = FindSparseIndexOf(index);
				CheckMutatedIndex(index, sparseIdx);
				values[sparseIdx] = value;
			}
		}

		/// <summary>
		/// Creates a new instance of <see cref="SparseVector"/> with the provided arrays as its internal data.
		/// </summary>
		/// <param name="length">The number of zero and non-zero entries of the new <see cref="SparseVector"/>.</param>
		/// <param name="values">The internal array that stores the values of the non-zero entries of the vector. Constraints: 
		///     <paramref name="values"/>.Length == <paramref name="indices"/>.Length &lt;= <paramref name="length"/>.</param>
		/// <param name="indices">The internal array that stores the indices of the non-zero entries of the vector. Constraints: 
		///     <paramref name="values"/>.Length == <paramref name="indices"/>.Length &lt;= <paramref name="length"/>.</param>
		/// <param name="checkInput">If true, the validity of <paramref name="values"/> and <paramref name="indices"/> will be 
		///     checked, which is safer. If false, no check will be made, which is daster.</param>
		/// <param name="sortInput">If true, <paramref name="values"/> and <paramref name="indices"/> will be sorted in
		///     ascending order of the entries of <paramref name="indices"/>. If false, they are assumed to be sorted. If they 
		///     are not, some methods may produce errors or have lower performance.</param>
		public static SparseVector CreateFromArrays(int length, double[] values, int[] indices,
			bool checkInput, bool sortInput)
		{
			bool verifiedSorted = false;
			if (sortInput)
			{
				Array.Sort(indices, values);
				verifiedSorted = true;
			}

			if (checkInput)
			{
				int nnz = indices.Length;
				if (values.Length != nnz)
				{
					throw new ArgumentException("The length of the values and indices arrays must be equal (and equal"
						+ $" to the number of non zero entries), but were {values.Length} and {indices.Length} respectively");
				}
				if (!verifiedSorted) // if they have already been sorted we do not need to check them again.
				{
					for (int i = 1; i < nnz; ++i)
					{
						if (indices[i] <= indices[i - 1]) throw new ArgumentException("The indices array must be sorted");
					}
				}
				int maxIndex = indices[nnz - 1]; // first sort them
				if (maxIndex >= length)
				{
					throw new ArgumentException($"This sparse vector contains {maxIndex + 1} entries, while its length is {length}");
				}
			}
			return new SparseVector(length, values, indices);
		}

		/// <summary>
		/// Creates a new instance of <see cref="SparseVector"/> with the entries of <paramref name="denseArray"/>. Only the
		/// entries that satisfy: <paramref name="denseArray"/>[i] != 0 are explicitly stored in the new 
		/// <see cref="SparseVector"/>.
		/// </summary>
		/// <param name="denseArray">The original vector that will be converted to <see cref="SparseVector"/>.</param>
		public static SparseVector CreateFromDense(double[] denseArray)
		{
			var indicesValues = new SortedDictionary<int, double>();
			for (int i = 0; i < denseArray.Length; ++i)
			{
				if (denseArray[i] != 0) indicesValues.Add(i, denseArray[i]);
			}
			return CreateFromDictionary(denseArray.Length, indicesValues);
		}

		/// <summary>
		/// Creates a new instance of <see cref="SparseVector"/> with the entries of <paramref name="denseArray"/>. Only the
		/// entries that satisfy: <paramref name="denseArray"/>[i] > <paramref name="tolerance"/> are explicitly stored in the 
		/// new <see cref="SparseVector"/>.
		/// </summary>
		/// <param name="denseArray">The original vector that will be converted to <see cref="SparseVector"/>.</param>
		/// <param name="tolerance">The tolerance under which an entry of <paramref name="denseArray"/> is considered to be 0. 
		///     Constraints: <paramref name="tolerance"/> &gt; 0. To keep only exact zeros, instead of setting 
		///     <paramref name="tolerance"/> = 0 use <see cref="CreateFromDense(double[])"/>.</param>
		public static SparseVector CreateFromDense(double[] denseArray, double tolerance)
		{
			var indicesValues = new SortedDictionary<int, double>();
			for (int i = 0; i < denseArray.Length; ++i)
			{
				if (Math.Abs(denseArray[i]) > tolerance) indicesValues.Add(i, denseArray[i]);
			}
			return CreateFromDictionary(denseArray.Length, indicesValues);
		}

		/// <summary>
		/// Creates a new instance of <see cref="SparseVector"/> with the entries of <paramref name="denseVector"/>. Only the
		/// entries that satisfy: <paramref name="denseVector"/>[i] != 0 are explicitly stored in the new 
		/// <see cref="SparseVector"/>.
		/// </summary>
		/// <param name="denseVector">The original vector that will be converted to <see cref="SparseVector"/>.</param>
		public static SparseVector CreateFromDense(Vector denseVector) => CreateFromDense(denseVector.RawData);

		/// <summary>
		/// Creates a new instance of <see cref="SparseVector"/> with the entries of <paramref name="denseVector"/>. Only the
		/// entries that satisfy: <paramref name="denseVector"/>[i] > <paramref name="tolerance"/> are explicitly stored in the 
		/// new <see cref="SparseVector"/>.
		/// </summary>
		/// <param name="denseVector">The original vector that will be converted to <see cref="SparseVector"/>.</param>
		/// <param name="tolerance">The tolerance under which an entry of <paramref name="denseVector"/> is considered to be 0. 
		///     Constraints: <paramref name="tolerance"/> &gt; 0. To keep only exact zeros, instead of setting 
		///     <paramref name="tolerance"/> = 0 use <see cref="CreateFromDense(double[])"/>.</param>
		public static SparseVector CreateFromDense(Vector denseVector, double tolerance)
		{
			return CreateFromDense(denseVector.RawData, tolerance);
		}

		/// <summary>
		/// Creates a new instance of <see cref="SparseVector"/> that has the provided <paramref name="length"/> and explicitly
		/// stores only the entries in <paramref name="nonZeroEntries"/>. All other entries are considered as 0. First 
		/// <paramref name="nonZeroEntries"/> will be sorted.
		/// </summary>
		/// <param name="length">The number of zero and non-zero entries in the new <see cref="SparseVector"/>.</param>
		/// <param name="nonZeroEntries">The indices and values of the non-zero entries of the new vector. Constraints:
		///     (foreach int idx in <paramref name="nonZeroEntries"/>.Keys) 0 &lt;= idx &lt; <paramref name="length"/>.</param>
		public static SparseVector CreateFromDictionary(int length, Dictionary<int, double> nonZeroEntries)
		{
			if (nonZeroEntries.Count == 0) new SparseVector(length, new double[0], new int[0]); // All zero vector

			double[] values = new double[nonZeroEntries.Count];
			int[] indices = new int[nonZeroEntries.Count];
			int nnz = 0;
			foreach (var idxValPair in nonZeroEntries.OrderBy(pair => pair.Key))
			{
				indices[nnz] = idxValPair.Key;
				values[nnz] = idxValPair.Value;
				++nnz;
			}

			int maxIndex = indices[nnz - 1]; // first sort them
			if (maxIndex >= length)
			{
				throw new ArgumentException($"This sparse vector contains {maxIndex + 1} entries, while its length is {length}");
			}
			return new SparseVector(length, values, indices);
		}

		/// <summary>
		/// Creates a new instance of <see cref="SparseVector"/> that has the provided <paramref name="length"/> and explicitly
		/// stores only the entries in <paramref name="nonZeroEntries"/>. All other entries are considered as 0.
		/// </summary>
		/// <param name="length">The number of zero and non-zero entries in the new <see cref="SparseVector"/>.</param>
		/// <param name="nonZeroEntries">The indices and values of the non-zero entries of the new vector. Constraints:
		///     (foreach int idx in <paramref name="nonZeroEntries"/>.Keys) 0 &lt;= idx &lt; <paramref name="length"/>.</param>
		public static SparseVector CreateFromDictionary(int length, SortedDictionary<int, double> nonZeroEntries)
		{
			if (nonZeroEntries.Count == 0) return new SparseVector(length, new double[0], new int[0]); // All zero vector

			double[] values = new double[nonZeroEntries.Count];
			int[] indices = new int[nonZeroEntries.Count];
			int nnz = 0;
			foreach (var idxValPair in nonZeroEntries)
			{
				indices[nnz] = idxValPair.Key;
				values[nnz] = idxValPair.Value;
				++nnz;
			}

			int maxIndex = indices[nnz - 1]; // first sort them
			if (maxIndex >= length)
			{
				throw new ArgumentException($"This sparse vector contains {maxIndex + 1} entries, while its length is {length}");
			}
			return new SparseVector(length, values, indices);
		}

		public override void AddToIndex(int index, double value)
		{
			int sparseIdx = FindSparseIndexOf(index);
			CheckMutatedIndex(index, sparseIdx);
			values[sparseIdx] = value;
		}

		public override IVector Axpy(IReadOnlyVector otherVector, double otherCoefficient)
		{
			Preconditions.CheckVectorDimensions(this, otherVector);
			if (otherVector is SparseVector otherSparse) // In case both matrices have the exact same index arrays
			{
				if (HasSameIndexer(otherSparse))
				{
					// Do not copy the index arrays, since they are already spread around. TODO: is this a good idea?
					double[] result = new double[this.values.Length];
					Array.Copy(this.values, result, this.values.Length);
					GlobalProvider.Blas.Daxpy(values.Length, otherCoefficient, otherSparse.values, 0, 1, result, 0, 1);
					return new SparseVector(Length, result, indices);
				}
			}
			else if (otherVector is Vector otherDense)
			{
				double[] result = otherDense.Scale(otherCoefficient).RawData;
				GlobalProvider.SparseBlas.Daxpyi(this.indices.Length, 1.0, this.values, this.indices, 0, result, 0);
				return Vector.CreateFromArray(result, false);
			}

			// All entries must be processed. TODO: optimizations may be possible (e.g. only access the nnz in this vector)
			return DenseStrategies.LinearCombination(this, 1.0, otherVector, otherCoefficient);
		}

		public override void AxpyIntoThis(IReadOnlyVector otherVector, double otherCoefficient)
		{
			if (otherVector is SparseVector otherSparse) AxpyIntoThis(otherSparse, otherCoefficient);
			else throw new SparsityPatternModifiedException(
				 "This operation is legal only if the other vector has the same sparsity pattern");
		}

		/// <summary>
		/// Performs the following operation for all i:
		/// this[i] = <paramref name="otherCoefficient"/> * <paramref name="otherVector"/>[i] + this[i]. 
		/// Optimized version of <see cref="IVector.DoEntrywise(IReadOnlyVector, Func{double, double, double})"/> and 
		/// <see cref="IVector.LinearCombination(double, IReadOnlyVector, double)"/>. Named after BLAS axpy (y = a*x plus y).
		/// The resulting vector overwrites the entries of this.
		/// </summary>
		/// <param name="otherVector">
		/// A vector with the same <see cref="IIndexable1D.Length"/> as this.
		/// </param>
		/// <param name="otherCoefficient">
		/// A scalar that multiplies each entry of <paramref name="otherVector"/>.
		/// </param>
		/// <exception cref="NonMatchingDimensionsException">
		/// Thrown if <paramref name="otherVector"/> has different <see cref="IIndexable1D.Length"/> than this.
		/// </exception>
		/// <exception cref="SparsityPatternModifiedException">
		/// Thrown if an entry this[i] needs to be overwritten, but that is not permitted by the vector storage format.
		/// </exception>
		public void AxpyIntoThis(SparseVector otherVector, double otherCoefficient)
		{
			Preconditions.CheckVectorDimensions(this, otherVector);
			if (!HasSameIndexer(otherVector)) throw new SparsityPatternModifiedException(
				"This operation is legal only if the other vector has the same sparsity pattern");

			GlobalProvider.Blas.Daxpy(values.Length, otherCoefficient, otherVector.values, 0, 1, this.values, 0, 1);
		}

		public override void AxpySubvectorIntoThis(int destinationIndex, IReadOnlyVector sourceVector, double sourceCoefficient,
			int sourceIndex, int length)
		{
			//TODO: needs testing for off-by-1 bugs and extension to cases where source and destination indices are different.

			Preconditions.CheckSubvectorDimensions(this, destinationIndex, length);
			Preconditions.CheckSubvectorDimensions(sourceVector, sourceIndex, length);

			if ((sourceVector is SparseVector otherSparse) && HasSameIndexer(otherSparse))
			{
				if (destinationIndex != sourceIndex) throw new NotImplementedException();
				int start = Array.FindIndex(this.indices, x => x >= destinationIndex);
				int end = Array.FindIndex(this.indices, x => x >= destinationIndex + length);
				int sparseLength = end - start;
				GlobalProvider.Blas.Daxpy(sparseLength, sourceCoefficient, otherSparse.values, start, 1, this.values, start, 1);
			}
			throw new SparsityPatternModifiedException(
				"This operation is legal only if the other vector has the same sparsity pattern");
		}

		public override void Clear() => Array.Clear(values, 0, values.Length);

		public override IVector Copy(bool copyIndexingData) => CopySparse();

		/// <summary>
		/// Initializes a new instance of <see cref="SparseVector"/> by deep copying the entries as this instance.
		/// </summary>
		public SparseVector CopySparse()
		{
			int n = values.Length;
			double[] valuesCopy = new double[n];
			Array.Copy(values, valuesCopy, n);
			int[] indicesCopy = new int[n];
			Array.Copy(indices, indicesCopy, n);
			return new SparseVector(n, valuesCopy, indicesCopy);
		}

		public override void CopyFrom(IReadOnlyVector sourceVector)
		{
			Preconditions.CheckVectorDimensions(this, sourceVector);
			if ((sourceVector is SparseVector otherSparse) && HasSameIndexer(otherSparse))
			{
				Array.Copy(otherSparse.values, this.values, this.Length);
			}

			throw new SparsityPatternModifiedException(
				 "This operation is legal only if the other vector has the same sparsity pattern");
		}

		public override void CopySubvectorFrom(int destinationIndex, IReadOnlyVector sourceVector, int sourceIndex, int length)
		{
			//TODO: needs testing for off-by-1 bugs and extension to cases where source and destination indices are different.
			Preconditions.CheckSubvectorDimensions(this, destinationIndex, length);
			Preconditions.CheckSubvectorDimensions(sourceVector, sourceIndex, length);

			if ((sourceVector is SparseVector otherSparse) && HasSameIndexer(otherSparse))
			{
				if (destinationIndex != sourceIndex) throw new NotImplementedException();
				int start = Array.FindIndex(this.indices, x => x >= destinationIndex);
				int end = Array.FindIndex(this.indices, x => x >= destinationIndex + length);
				int sparseLength = end - start;
				Array.Copy(otherSparse.values, start, this.values, start, sparseLength);
			}
			throw new SparsityPatternModifiedException(
				"This operation is legal only if the other vector has the same sparsity pattern");
		}

		public override double[] CopyToArray()
		{
			double[] result = new double[Length];
			for (int i = 0; i < values.Length; ++i) result[indices[i]] = values[i];
			return result;
		}

		/// <summary>
		/// Initializes a new instance of <see cref="Vector"/> that contains the same non-zero entries as this 
		/// <see cref="SparseVector"/>, while the rest entries are explicitly 0.
		/// </summary>
		public Vector CopyToFullVector() => Vector.CreateFromArray(CopyToArray(), false);

		/// <summary>
		/// Counts how many non zero entries are stored in the vector. This includes zeros that are explicitly stored.
		/// </summary>
		public int CountNonZeros() => values.Length;

		public override IVector CreateZeroVectorWithSameFormat() => CreateZeroVector();

		public SparseVector CreateZeroVector() => new SparseVector(Length, new double[indices.Length], indices);

		public override IVector DoEntrywise(IReadOnlyVector otherVector, Func<double, double, double> binaryOperation)
		{
			Preconditions.CheckVectorDimensions(this, otherVector);
			if (otherVector is SparseVector otherSparse) // In case both matrices have the exact same index arrays
			{
				if (HasSameIndexer(otherSparse))
				{
					// Do not copy the index arrays, since they are already spread around. TODO: is this a good idea?
					double[] resultValues = new double[values.Length];
					for (int i = 0; i < values.Length; ++i)
					{
						resultValues[i] = binaryOperation(this.values[i], otherSparse.values[i]);
					}
					return new SparseVector(Length, resultValues, indices);
				}
			}

			// All entries must be processed. TODO: optimizations may be possible (e.g. only access the nnz in this vector)
			return DenseStrategies.DoEntrywise(this, otherVector, binaryOperation);
		}

		public override void DoEntrywiseIntoThis(IReadOnlyVector otherVector, Func<double, double, double> binaryOperation)
		{
			Preconditions.CheckVectorDimensions(this, otherVector);
			if ((otherVector is SparseVector otherSparse) && HasSameIndexer(otherSparse))
			{
				for (int i = 0; i < values.Length; ++i)
				{
					this.values[i] = binaryOperation(this.values[i], otherSparse.values[i]);
				}
			}
			throw new SparsityPatternModifiedException(
				 "This operation is legal only if the other vector has the same sparsity pattern");
		}

		public override IVector DoToAllEntries(Func<double, double> unaryOperation)
		{
			// Only apply the operation on non zero entries
			double[] newValues = new double[values.Length];
			for (int i = 0; i < values.Length; ++i) newValues[i] = unaryOperation(values[i]);

			if (new ValueComparer(1e-10).AreEqual(unaryOperation(0.0), 0.0)) // The same sparsity pattern can be used.
			{
				// Copy the index arrays. TODO: See if we can use the same index arrays (e.g. if this class does not change them (it shouldn't))
				int[] indicesCopy = new int[indices.Length];
				Array.Copy(indices, indicesCopy, indices.Length);
				return new SparseVector(Length, newValues, indicesCopy);
			}
			else // The sparsity is destroyed. Revert to a full vector.
			{
				return new SparseVector(Length, newValues, indices).CopyToFullVector();
			}
		}

		public override void DoToAllEntriesIntoThis(Func<double, double> unaryOperation)
		{
			if (new ValueComparer(1e-10).AreEqual(unaryOperation(0.0), 0.0))
			{
				for (int i = 0; i < values.Length; ++i) values[i] = unaryOperation(values[i]);
			}
			else throw new SparsityPatternModifiedException("This operation will change the sparsity pattern");
		}

		public override double DotProduct(IReadOnlyVector vector)
		{
			Preconditions.CheckVectorDimensions(this, vector);

			if (vector is Vector dense)
			{
				return GlobalProvider.SparseBlas.Ddoti(values.Length, values, indices, 0, dense.RawData, 0);
			}
			else if ((vector is SparseVector sparse) && HasSameIndexer(sparse))
			{
				return GlobalProvider.Blas.Ddot(values.Length, this.values, 0, 1, sparse.values, 0, 1);
			}

			double sum = 0;
			for (int i = 0; i < values.Length; ++i) sum += values[i] * vector[indices[i]];
			return sum;
		}

		public override bool Equals(IIndexable1D other, double tolerance = 1e-13)
		{
			if (this.Length != other.Length) return false;
			var comparer = new ValueComparer(tolerance);
			int previousIndex = 0;
			for (int i = 0; i < indices.Length; ++i)
			{
				int index = indices[i];
				for (int j = previousIndex; j < index; ++j) // zero entries between the stored ones
				{
					if (!comparer.AreEqual(0.0, other[j])) return false;
				}
				if (!comparer.AreEqual(values[i], other[index])) return false; // Non zero entry
				previousIndex = index + 1;
			}
			return true;
		}

		/// <summary>
		/// Iterates over the non zero entries of the vector. This includes zeros that are explicitly stored.
		/// </summary>
		public IEnumerable<(int index, double value)> EnumerateNonZeros()
		{
			for (int i = 0; i < values.Length; ++i)
			{
				yield return (indices[i], values[i]);
			}
		}

		public override bool HasSameFormat(IReadOnlyVector other)
		{
			if (other is SparseVector casted && casted.indices == this.indices)
			{
				return true;
			}
			else
			{
				return false;
			}
		}

		public override IVector LinearCombination(double thisCoefficient, IReadOnlyVector otherVector, double otherCoefficient)
		{
			Preconditions.CheckVectorDimensions(this, otherVector);
			if (otherVector is SparseVector otherSparse) // In case both matrices have the exact same index arrays
			{
				if (HasSameIndexer(otherSparse))
				{
					// Do not copy the index arrays, since they are already spread around. TODO: is this a good idea?
					double[] result = new double[this.values.Length];
					if (thisCoefficient == 1.0)
					{
						Array.Copy(this.values, result, this.values.Length);
						GlobalProvider.Blas.Daxpy(values.Length, otherCoefficient, otherSparse.values, 0, 1, result, 0, 1);
					}
					else if (otherCoefficient == 1.0)
					{
						Array.Copy(otherSparse.values, result, values.Length);
						GlobalProvider.Blas.Daxpy(values.Length, thisCoefficient, this.values, 0, 1, result, 0, 1);
					}
					else
					{
						Array.Copy(this.values, result, this.values.Length);
						GlobalProvider.Blas.Daxpby(values.Length, otherCoefficient, otherSparse.values, 0, 1,
							thisCoefficient, result, 0, 1);
					}
					return new SparseVector(Length, result, indices);
				}
			}

			// All entries must be processed. TODO: optimizations may be possible (e.g. only access the nnz in this vector)
			return DenseStrategies.LinearCombination(this, thisCoefficient, otherVector, otherCoefficient);
		}

		public override void LinearCombinationIntoThis(double thisCoefficient, IReadOnlyVector otherVector, double otherCoefficient)
		{
			Preconditions.CheckVectorDimensions(this, otherVector);
			if ((otherVector is SparseVector otherSparse) && HasSameIndexer(otherSparse))
			{
				if (thisCoefficient == 1.0)
				{
					GlobalProvider.Blas.Daxpy(values.Length, otherCoefficient, otherSparse.values, 0, 1, this.values, 0, 1);
				}
				else
				{
					GlobalProvider.Blas.Daxpby(values.Length, otherCoefficient, otherSparse.values, 0, 1,
						thisCoefficient, this.values, 0, 1);
				}
			}

			throw new SparsityPatternModifiedException(
				 "This operation is legal only if the other vector has the same sparsity pattern");
		}

		public override double Norm2() => GlobalProvider.Blas.Dnrm2(values.Length, values, 0, 1);

		public override double Reduce(
			double identityValue, ProcessEntry processEntry, ProcessZeros processZeros, Finalize finalize)
		{
			double aggregator = identityValue;
			int nnz = values.Length;
			for (int i = 0; i < nnz; ++i) aggregator = processEntry(values[i], aggregator);
			aggregator = processZeros(Length - nnz, aggregator);
			return finalize(aggregator);
		}

		public override IVector Scale(double scalar) => ScaleSparse(scalar);

		/// <summary>
		/// Performs the following operation for 0 &lt;= i &lt; this.<see cref="Length"/>:
		/// result[i] = <paramref name="scalar"/> * this[i].
		/// The resulting vector is written to a new <see cref="SparseVector"/> and then returned.
		/// </summary>
		/// <param name="scalar">The scalar value that multiplies all entries of the vector.</param>
		public SparseVector ScaleSparse(double scalar)
		{
			int nnz = this.values.Length;
			double[] resultValues = new double[nnz];
			Array.Copy(this.values, resultValues, nnz);
			GlobalProvider.Blas.Dscal(nnz, scalar, resultValues, 0, 1);
			return new SparseVector(Length, resultValues, this.indices); //TODO: perhaps I should also copy the indices
		}

		public override void ScaleIntoThis(double scalar) => GlobalProvider.Blas.Dscal(values.Length, scalar, values, 0, 1);

		public override void SetAll(double value)
		{
			if (value == 0.0)
			{
				Clear();
			}
			else
			{
				throw new SparsityPatternModifiedException("Cannot change structural zero entries.");
			}
		}

		/// <summary>
		/// Returns the index into <see cref="values"/> of the entry this[<paramref name="denseIdx"/>]. If this entry is a
		/// structural zero, -1 will be returned.
		/// </summary>
		/// <param name="denseIdx"></param>
		private int FindSparseIndexOf(int denseIdx)
		{
			Preconditions.CheckIndex1D(this, denseIdx);
			return Array.BinarySearch(indices, denseIdx); // only works if indices are sorted!!!
		}

		private bool HasSameIndexer(SparseVector other)
		{
			return this.indices == other.indices;
		}

		private void CheckMutatedIndex(int index, int sparseIdx)
		{
			if (sparseIdx < 0) throw new SparsityPatternModifiedException(
				$"The entry at index = {index} is zero and not stored explicilty, therefore it cannot be modified.");
		}
	}
}
