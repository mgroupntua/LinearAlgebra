namespace MGroup.LinearAlgebra.Vectors
{
	using System;
	using System.Diagnostics;

	using MGroup.LinearAlgebra.Commons;
	using MGroup.LinearAlgebra.Exceptions;
	using MGroup.LinearAlgebra.Reduction;

	public abstract class DefaultVector : IVector
	{
		public abstract int Length { get; }

		public abstract double this[int index] { get; }

		public abstract void Clear();

		public abstract IVector CreateZeroVectorWithSameFormat();

		public abstract bool HasSameFormat(IVectorView other);

		public abstract void Set(int index, double value);

		public virtual void AddIntoThisNonContiguouslyFrom(int[] thisIndices, IVectorView otherVector, int[] otherIndices)
		{
			if (thisIndices.Length != otherIndices.Length)
			{
				throw new NonMatchingDimensionsException("Must operate on the same number of indices on both vectors");
			}

			WarnAboutPerformanceBottlenecks();
			ProhibitPerformanceBottlenecks();
			for (int i = 0; i < otherIndices.Length; ++i)
			{
				double val = otherVector[otherIndices[i]];
				this.Set(thisIndices[i], this[i] + val);
			}
		}

		public virtual void AddIntoThisNonContiguouslyFrom(int[] thisIndices, IVectorView otherVector)
		{
			WarnAboutPerformanceBottlenecks();
			ProhibitPerformanceBottlenecks();
			for (int i = 0; i < thisIndices.Length; ++i)
			{
				this.Set(thisIndices[i], otherVector[i]);
			}
		}

		public virtual void AddToIndex(int index, double value) => Set(index, this[index] + value);

		public virtual IVector Axpy(IVectorView otherVector, double otherCoefficient)
			=> LinearCombination(1.0, otherVector, otherCoefficient);

		public virtual void AxpyIntoThis(IVectorView otherVector, double otherCoefficient)
			=> LinearCombinationIntoThis(1.0, otherVector, otherCoefficient);

		public virtual void AxpySubvectorIntoThis(
			int destinationIndex, IVectorView sourceVector, double sourceCoefficient, int sourceIndex, int length)
		{
			if (destinationIndex + length <= this.Length)
			{
				throw new NonMatchingDimensionsException("Not enough space on this vector to write the requested entries.");
			}

			if (sourceIndex + length <= sourceVector.Length)
			{
				throw new NonMatchingDimensionsException("The source vector does not have as many entries as requested.");
			}

			WarnAboutPerformanceBottlenecks();
			ProhibitPerformanceBottlenecks();
			for (int i = 0; i < length; ++i)
			{
				double val = sourceCoefficient * sourceVector[sourceIndex + i];
				this.Set(destinationIndex + i, this[destinationIndex + i] + val);
			}
		}

		public virtual IVector Copy(bool copyIndexingData = false)
		{
			WarnAboutPerformanceBottlenecks();
			ProhibitPerformanceBottlenecks();
			IVector clone = CreateZeroVectorWithSameFormat();
			for (int i = 0; i < Length; i++)
			{
				double val = this[i];
				if (val != 0.0)
				{
					clone.Set(i, val);
				}
			}

			return clone;
		}

		public virtual void CopyFrom(IVectorView sourceVector)
		{
			Preconditions.CheckVectorDimensions(this, sourceVector);
			WarnAboutPerformanceBottlenecks();
			ProhibitPerformanceBottlenecks();
			Clear();
			for (int i = 0; i < Length; i++)
			{
				double val = sourceVector[i];
				if (val != 0.0)
				{
					this.Set(i, val);
				}
			}
		}

		public virtual void CopyNonContiguouslyFrom(int[] thisIndices, IVectorView otherVector, int[] otherIndices)
		{
			if (thisIndices.Length != otherIndices.Length)
			{
				throw new NonMatchingDimensionsException("Must operate on the same number of indices on both vectors");
			}

			WarnAboutPerformanceBottlenecks();
			ProhibitPerformanceBottlenecks();
			for (int i = 0; i < otherIndices.Length; ++i)
			{
				double val = otherVector[otherIndices[i]];
				this.Set(thisIndices[i], val);
			}
		}

		public virtual void CopyNonContiguouslyFrom(IVectorView otherVector, int[] otherIndices)
		{
			WarnAboutPerformanceBottlenecks();
			ProhibitPerformanceBottlenecks();
			for (int i = 0; i < otherIndices.Length; ++i)
			{
				double val = otherVector[otherIndices[i]];
				this.Set(i, val);
			}
		}

		public virtual void CopySubvectorFrom(int destinationIndex, IVectorView sourceVector, int sourceIndex, int length)
		{
			if (destinationIndex + length <= this.Length)
			{
				throw new NonMatchingDimensionsException("Not enough space on this vector to write the requested entries.");
			}

			if (sourceIndex + length <= sourceVector.Length)
			{
				throw new NonMatchingDimensionsException("The source vector does not have as many entries as requested.");
			}

			WarnAboutPerformanceBottlenecks();
			ProhibitPerformanceBottlenecks();
			for (int i = 0; i < length; ++i)
			{
				double val = sourceVector[sourceIndex + i];
				this.Set(destinationIndex + i, val);
			}
		}

		public virtual double[] CopyToArray()
		{
			WarnAboutPerformanceBottlenecks();
			ProhibitPerformanceBottlenecks();
			var result = new double[Length];
			for (int i = 0; i < result.Length; i++)
			{
				result[i] = this[i];
			}

			return result;
		}

		public virtual IVector DoEntrywise(IVectorView otherVector, Func<double, double, double> binaryOperation)
		{
			Preconditions.CheckVectorDimensions(this, otherVector);
			if (this.HasSameFormat(otherVector) && (binaryOperation(0.0, 0.0) == 0.0))
			{
				IVector result = Copy(copyIndexingData: false);
				result.DoEntrywiseIntoThis(otherVector, binaryOperation);
				return result;
			}
			else
			{
				WarnAboutPerformanceBottlenecks();
				ProhibitPerformanceBottlenecks();
				var result = new double[Length];
				for (int i = 0; i < Length; i++)
				{
					result[i] = binaryOperation(this[i], otherVector[i]);
				}

				return Vector.CreateFromArray(result, false);
			}
		}

		public virtual void DoEntrywiseIntoThis(IVectorView otherVector, Func<double, double, double> binaryOperation)
		{
			Preconditions.CheckVectorDimensions(this, otherVector);
			WarnAboutPerformanceBottlenecks();
			ProhibitPerformanceBottlenecks();
			for (int i = 0; i < Length; i++)
			{
				this.Set(i, binaryOperation(this[i], otherVector[i]));
			}
		}

		public virtual IVector DoToAllEntries(Func<double, double> unaryOperation)
		{
			WarnAboutPerformanceBottlenecks();
			ProhibitPerformanceBottlenecks();
			IVector result;
			if (unaryOperation(0.0) == 0.0)
			{
				result = Copy(copyIndexingData: false);
			}
			else
			{
				result = Vector.CreateFromArray(this.CopyToArray(), false);
			}

			result.DoToAllEntriesIntoThis(unaryOperation);
			return result;
		}

		public virtual void DoToAllEntriesIntoThis(Func<double, double> unaryOperation)
		{
			WarnAboutPerformanceBottlenecks();
			ProhibitPerformanceBottlenecks();
			for (int i = 0; i < Length; i++)
			{
				this.Set(i, unaryOperation(this[i]));
			}
		}

		public virtual double DotProduct(IVectorView vector)
		{
			Preconditions.CheckVectorDimensions(this, vector);
			WarnAboutPerformanceBottlenecks();
			ProhibitPerformanceBottlenecks();
			double sum = 0;
			for (int i = 0; i < Length; i++)
			{
				sum += this[i] * vector[i];
			}

			return sum;
		}

		public virtual bool Equals(IIndexable1D other, double tolerance = 1E-13)
		{
			if (this.Length != other.Length)
			{
				return false;
			}

			WarnAboutPerformanceBottlenecks();
			ProhibitPerformanceBottlenecks();

			var comparer = new ValueComparer(tolerance);
			for (int i = 0; i < Length; ++i)
			{
				if (!comparer.AreEqual(this[i], other[i]))
				{
					return false;
				}
			}

			return true;
		}

		public virtual IVector LinearCombination(double thisCoefficient, IVectorView otherVector, double otherCoefficient)
			=> DoEntrywise(otherVector, (x, y) => thisCoefficient * x + otherCoefficient * y);

		public virtual void LinearCombinationIntoThis(double thisCoefficient, IVectorView otherVector, double otherCoefficient)
			=> DoEntrywiseIntoThis(otherVector, (x, y) => thisCoefficient * x + otherCoefficient * y);

		public virtual double Norm2() => Reductions.Norm2(this);

		public virtual double Reduce(double identityValue, ProcessEntry processEntry, ProcessZeros processZeros, Finalize finalize)
		{
			WarnAboutPerformanceBottlenecks();
			ProhibitPerformanceBottlenecks();
			double accumulator = identityValue;
			for (int i = 0; i < this.Length; ++i)
			{
				accumulator = processEntry(this[i], accumulator);
			}

			return finalize(accumulator);
		}

		public virtual IVector Scale(double scalar) => DoToAllEntries(x => scalar * x);

		public virtual void ScaleIntoThis(double scalar) => DoToAllEntries(x => scalar * x);

		public virtual void SetAll(double value)
		{
			WarnAboutPerformanceBottlenecks();
			ProhibitPerformanceBottlenecks();
			for (int i = 0; i < Length; i++)
			{
				Set(i, value);
			}
		}

		[Conditional("DEBUG")]
		private static void WarnAboutPerformanceBottlenecks()
		{
			Debug.WriteLine(
				"Potential performance bottleneck due to accessing all entries of a potentially sparse matrix or vector, at :"
				+ Environment.StackTrace);
		}

		[Conditional("RELEASE")]
		private static void ProhibitPerformanceBottlenecks()
		{
			if (LibrarySettings.ThrowExceptionOnKnownPerformanceBottlenecksInReleaseBuilds)
			{
				throw new PerformanceBottleneckException(
					"Potential performance bottleneck due to accessing all entries of a potentially sparse matrix or vector.");
			}
		}
	}
}
