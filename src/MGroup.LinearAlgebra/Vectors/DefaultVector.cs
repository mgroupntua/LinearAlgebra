namespace MGroup.LinearAlgebra.Vectors
{
	using System;
	using System.Diagnostics;

	using MGroup.LinearAlgebra.Commons;
	using MGroup.LinearAlgebra.Exceptions;
	using MGroup.LinearAlgebra.Reduction;

	using static MGroup.LinearAlgebra.Commons.PerformanceWarnings;

	public abstract class DefaultVector : IVector
	{
		public abstract int Length { get; }

		public abstract double this[int index] { get; set; }

		public abstract void Clear();

		public abstract IVector CreateZeroVectorWithSameFormat();

		public abstract bool HasSameFormat(IReadOnlyVector other);

		public virtual void AddIntoThisNonContiguouslyFrom(int[] thisIndices, IReadOnlyVector otherVector, int[] otherIndices)
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
				this[thisIndices[i]] += val;
			}
		}

		public virtual void AddIntoThisNonContiguouslyFrom(int[] thisIndices, IReadOnlyVector otherVector)
		{
			WarnAboutPerformanceBottlenecks();
			ProhibitPerformanceBottlenecks();
			for (int i = 0; i < thisIndices.Length; ++i)
			{
				this[thisIndices[i]] += otherVector[i];
			}
		}

		public virtual void AddToIndex(int index, double value) => this[index] = this[index] + value;

		public virtual IVector Axpy(IReadOnlyVector otherVector, double otherCoefficient)
			=> LinearCombination(1.0, otherVector, otherCoefficient);

		public virtual void AxpyIntoThis(IReadOnlyVector otherVector, double otherCoefficient)
			=> LinearCombinationIntoThis(1.0, otherVector, otherCoefficient);

		public virtual void AxpySubvectorIntoThis(
			int destinationIndex, IReadOnlyVector sourceVector, double sourceCoefficient, int sourceIndex, int length)
		{
			if (destinationIndex + length > this.Length)
			{
				throw new NonMatchingDimensionsException("Not enough space on this vector to write the requested entries.");
			}

			if (sourceIndex + length > sourceVector.Length)
			{
				throw new NonMatchingDimensionsException("The source vector does not have as many entries as requested.");
			}

			WarnAboutPerformanceBottlenecks();
			ProhibitPerformanceBottlenecks();
			for (int i = 0; i < length; ++i)
			{
				this[destinationIndex + i] += sourceCoefficient * sourceVector[sourceIndex + i];
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

		public virtual void CopyFrom(IReadOnlyVector sourceVector)
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
					this[i] = val;
				}
			}
		}

		public virtual void CopyNonContiguouslyFrom(int[] thisIndices, IReadOnlyVector otherVector, int[] otherIndices)
		{
			if (thisIndices.Length != otherIndices.Length)
			{
				throw new NonMatchingDimensionsException("Must operate on the same number of indices on both vectors");
			}

			WarnAboutPerformanceBottlenecks();
			ProhibitPerformanceBottlenecks();
			for (int i = 0; i < otherIndices.Length; ++i)
			{
				this[thisIndices[i]] = otherVector[otherIndices[i]];
			}
		}

		public virtual void CopyNonContiguouslyFrom(IReadOnlyVector otherVector, int[] otherIndices)
		{
			WarnAboutPerformanceBottlenecks();
			ProhibitPerformanceBottlenecks();
			for (int i = 0; i < otherIndices.Length; ++i)
			{
				this[i] = otherVector[otherIndices[i]];
			}
		}

		public virtual void CopySubvectorFrom(int destinationIndex, IReadOnlyVector sourceVector, int sourceIndex, int length)
		{
			if (destinationIndex + length > this.Length)
			{
				throw new NonMatchingDimensionsException("Not enough space on this vector to write the requested entries.");
			}

			if (sourceIndex + length > sourceVector.Length)
			{
				throw new NonMatchingDimensionsException("The source vector does not have as many entries as requested.");
			}

			WarnAboutPerformanceBottlenecks();
			ProhibitPerformanceBottlenecks();
			for (int i = 0; i < length; ++i)
			{
				this[destinationIndex + i] = sourceVector[sourceIndex + i];
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

		public virtual IVector DoEntrywise(IReadOnlyVector otherVector, Func<double, double, double> binaryOperation)
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

		public virtual void DoEntrywiseIntoThis(IReadOnlyVector otherVector, Func<double, double, double> binaryOperation)
		{
			Preconditions.CheckVectorDimensions(this, otherVector);
			WarnAboutPerformanceBottlenecks();
			ProhibitPerformanceBottlenecks();
			for (int i = 0; i < Length; i++)
			{
				this[i] = binaryOperation(this[i], otherVector[i]);
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
				this[i] = unaryOperation(this[i]);
			}
		}

		public virtual double DotProduct(IReadOnlyVector vector)
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

		public virtual IVector LinearCombination(double thisCoefficient, IReadOnlyVector otherVector, double otherCoefficient)
		{
			Preconditions.CheckVectorDimensions(this, otherVector);
			if (this.HasSameFormat(otherVector))
			{
				IVector result = Copy(copyIndexingData: false);
				result.LinearCombinationIntoThis(thisCoefficient, otherVector, otherCoefficient);
				return result;
			}
			else
			{
				WarnAboutPerformanceBottlenecks();
				ProhibitPerformanceBottlenecks();
				var result = new double[Length];
				for (int i = 0; i < Length; i++)
				{
					result[i] = thisCoefficient * this[i] + otherCoefficient * otherVector[i];
				}

				return Vector.CreateFromArray(result, false);
			}
		}

		public virtual void LinearCombinationIntoThis(double thisCoefficient, IReadOnlyVector otherVector, double otherCoefficient)
		{
			Preconditions.CheckVectorDimensions(this, otherVector);
			WarnAboutPerformanceBottlenecks();
			ProhibitPerformanceBottlenecks();
			for (int i = 0; i < Length; i++)
			{
				this[i] = thisCoefficient * this[i] + otherCoefficient * otherVector[i];
			}
		}

		public virtual double Norm2()
		{
			WarnAboutPerformanceBottlenecks();
			ProhibitPerformanceBottlenecks();
			double sum = 0.0;
			for (int i = 0; i < this.Length; ++i)
			{
				double x = this[i];
				sum += x * x;
			}

			return Math.Sqrt(sum);
		}

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

		public virtual IVector Scale(double scalar)
		{
			WarnAboutPerformanceBottlenecks();
			ProhibitPerformanceBottlenecks();
			IVector result;
			result = Copy(copyIndexingData: false);
			result.ScaleIntoThis(scalar);
			return result;
		}

		public virtual void ScaleIntoThis(double scalar)
		{
			WarnAboutPerformanceBottlenecks();
			ProhibitPerformanceBottlenecks();
			for (int i = 0; i < Length; i++)
			{
				this[i] *= scalar;
			}
		}

		public virtual void Set(int index, double value) => this[index] = value;

		public virtual void SetAll(double value)
		{
			WarnAboutPerformanceBottlenecks();
			ProhibitPerformanceBottlenecks();
			for (int i = 0; i < Length; i++)
			{
				Set(i, value);
			}
		}
	}
}
