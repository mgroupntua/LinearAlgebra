using System;
using MGroup.LinearAlgebra.Commons;
using MGroup.LinearAlgebra.Exceptions;
using MGroup.LinearAlgebra.Reduction;

namespace MGroup.LinearAlgebra.Vectors
{
    /// <summary>
    /// A vector with 2 entries. Optimized version of <see cref="Vector"/>.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class Vector2: IVector, IEntrywiseOperableView1D<Vector2, Vector2>, IEntrywiseOperable1D<Vector2>
    {
        private readonly double[] data;

        private Vector2(double[] data)
        {
            this.data = data;
        }

        /// <summary>
        /// See <see cref="IIndexable1D.Length"/>.
        /// </summary>
        public int Length { get { return 2; } }

        /// <summary>
        /// The internal array that stores the entries of the vector. It should only be used for passing the raw array to linear 
        /// algebra libraries.
        /// </summary>
        internal double[] InternalData { get { return data; } }

        /// <summary>
        /// See <see cref="IIndexable1D.this[int]"/>.
        /// </summary>
        public double this[int index]
        {
            get { return data[index]; }
            set { data[index] = value; }
        }

        /// <summary>
        /// Initializes a new instance of <see cref="Vector2"/> that contains the provided entries: 
        /// {<paramref name="entry0"/>, <paramref name="entry1"/>}.
        /// </summary>
        /// <param name="entry0">The first entry of the new vector.</param>
        /// <param name="entry1">The second entry of the new vector.</param>
        public static Vector2 Create(double entry0, double entry1) => new Vector2(new double[] { entry0, entry1 });

        /// <summary>
        /// Initializes a new instance of <see cref="Vector2"/> with <paramref name="data"/> or a clone as its internal array.
        /// </summary>
        /// <param name="data">The entries of the vector to create. Constraints: <paramref name="data"/>.Length == 2.</param>
        /// <param name="copyArray">If true, <paramref name="data"/> will be copied and the new <see cref="Vector2"/> instance 
        ///     will have a reference to the copy, which is safer. If false, the new vector will have a reference to 
        ///     <paramref name="data"/> itself, which is faster.</param>
        /// <returns></returns>
        public static Vector2 CreateFromArray(double[] data, bool copyArray = false)
        {
            if (data.Length != 2) throw new NonMatchingDimensionsException(
                $"The provided array had length = {data.Length} instead of 2");
            if (copyArray) return new Vector2(new double[] { data[0], data[1] });
            else return new Vector2(data);
        }

        /// <summary>
        /// Creates an new instance of <see cref="Vector2"/> with all entries being equal to 0.
        /// </summary>
        public static Vector2 CreateZero() => new Vector2(new double[2]);

        #region operators 
        /// <summary>
        /// Performs the operation: result[i] = <paramref name="v1"/>[i] + <paramref name="v2"/>[i], 
        /// for 0 &lt;= i &lt; 2. The resulting entries are written to a new <see cref="Vector2"/> instance.
        /// </summary>
        /// <param name="v1">The first <see cref="Vector2"/> operand.</param>
        /// <param name="v2">The second <see cref="Vector2"/> operand.</param>
        public static Vector2 operator +(Vector2 v1, Vector2 v2)
        {
            return new Vector2(new double[] {v1.data[0] + v2.data[0], v1.data[1] + v2.data[1] });
        }

        /// <summary>
        /// Performs the operation: result[i] = <paramref name="v1"/>[i] - <paramref name="v2"/>[i], 
        /// for 0 &lt;= i &lt; 2. The resulting entries are written to a new <see cref="Vector2"/> instance.
        /// </summary>
        /// <param name="v1">The first <see cref="Vector2"/> operand.</param>
        /// <param name="v2">The second <see cref="Vector2"/> operand.</param>
        public static Vector2 operator -(Vector2 v1, Vector2 v2)
        {
            return new Vector2(new double[] { v1.data[0] - v2.data[0], v1.data[1] - v2.data[1] });
        }

        /// <summary>
        /// Performs the operation: result[i] = <paramref name="scalar"/> * <paramref name="vector"/>[i],
        /// for 0 &lt;= i &lt; 2. The resulting entries are written to a new <see cref="Vector2"/> instance.
        /// </summary>
        /// <param name="scalar">The scalar value that will be multiplied with all vector entries.</param>
        /// <param name="vector">The vector to multiply.</param>
        public static Vector2 operator *(double scalar, Vector2 vector)
        {
            return new Vector2(new double[] { scalar * vector.data[0], scalar * vector.data[1] });
        }

        /// <summary>
        /// Performs the operation: result[i] = <paramref name="scalar"/> * <paramref name="vector"/>[i],
        /// for 0 &lt;= i &lt; 2. The resulting entries are written to a new <see cref="Vector2"/> instance.
        /// </summary>
        /// <param name="vector">The vector to multiply.</param>
        /// <param name="scalar">The scalar value that will be multiplied with all vector entries.</param>
        public static Vector2 operator *(Vector2 vector, double scalar)
        {
            return new Vector2(new double[] { scalar * vector.data[0], scalar * vector.data[1] });
        }

        /// <summary>
        /// Performs the operation: scalar = sum(<paramref name="v1"/>[i] * <paramref name="v2"/>[i]), 
        /// for 0 &lt;= i &lt; 2. The resulting entries are written to a new <see cref="Vector2"/> instance.
        /// </summary>
        /// <param name="v1">The first <see cref="Vector2"/> operand.</param>
        /// <param name="v2">The second <see cref="Vector2"/> operand.</param>
        public static double operator *(Vector2 v1, Vector2 v2)
        {
            return v1.data[0] * v2.data[0] + v1.data[1] * v2.data[1];
        }
        #endregion

        /// <summary>
        /// Performs the operation: this[i] = this[i] + <paramref name="vector"/>[i], 
        /// for 0 &lt;= i &lt; 2. The resulting vector overwrites the entries of this <see cref="Vector2"/> instance.
        /// </summary>
        /// <param name="vector">A vector with two entries.</param>
        public void AddIntoThis(Vector2 vector)
        {
            this.data[0] += vector.data[0];
            this.data[1] += vector.data[1];
        }

        /// <summary>
        /// See <see cref="IVector.AddIntoThisNonContiguouslyFrom(int[], IVectorView, int[])"/>
        /// </summary>
        public void AddIntoThisNonContiguouslyFrom(int[] thisIndices, IVectorView otherVector, int[] otherIndices)
            => DenseStrategies.AddNonContiguouslyFrom(this, thisIndices, otherVector, otherIndices);

        /// <summary>
        /// See <see cref="IVector.AddIntoThisNonContiguouslyFrom(int[], IVectorView)"/>
        /// </summary>
        public void AddIntoThisNonContiguouslyFrom(int[] thisIndices, IVectorView otherVector)
            => DenseStrategies.AddNonContiguouslyFrom(this, thisIndices, otherVector);

		/// <summary>
		/// See <see cref="IVector.AddToIndex(int, double)"/>
		/// </summary>
		public void AddToIndex(int index, double value)
		{
			data[index] += value;
		}

		/// <summary>
		/// See <see cref="IVectorView.Axpy(IVectorView, double)"/>.
		/// </summary>
		public IVector Axpy(IVectorView otherVector, double otherCoefficient)
        {
            if (otherVector is Vector2 casted) return Axpy(casted, otherCoefficient);
            else
            {
                Preconditions.CheckVectorDimensions(this, otherVector);
                return new Vector2(new double[]
                {
                    data[0] + otherCoefficient * otherVector[0], data[1] + otherCoefficient * otherVector[1]
                });
            }
        }

        /// <summary>
        /// Performs the operation: result[i] = <paramref name="otherCoefficient"/> * <paramref name="otherVector"/>[i] + this[i], 
        /// for 0 &lt;= i &lt; 2. The resulting vector is written to a new <see cref="Vector2"/> and then returned.
        /// </summary>
        /// <param name="otherVector">A vector with two entries.</param>
        /// <param name="otherCoefficient">A scalar that multiplies each entry of <paramref name="otherVector"/>.</param>
        public Vector2 Axpy(Vector2 otherVector, double otherCoefficient)
        {
            return new Vector2(new double[] 
            {
                this.data[0] + otherCoefficient * otherVector.data[0], this.data[1] + otherCoefficient * otherVector.data[1]
            });
        }

        /// <summary>
        /// See <see cref="IVector.AxpyIntoThis(IVectorView, double)"/>
        /// </summary>
        public void AxpyIntoThis(IVectorView otherVector, double otherCoefficient)
        {
            if (otherVector is Vector2 casted) AxpyIntoThis(casted, otherCoefficient);
            else
            {
                Preconditions.CheckVectorDimensions(this, otherVector);
                data[0] += otherCoefficient * otherVector[0];
                data[1] += otherCoefficient * otherVector[1];
            }
        }

        /// <summary>
        /// Performs the operation: this[i] = <paramref name="otherCoefficient"/> * <paramref name="otherVector"/>[i] + this[i], 
        /// for 0 &lt;= i &lt; 2. The resulting vector overwrites the entries of this <see cref="Vector2"/> instance.
        /// </summary>
        /// <param name="otherVector">A vector with two entries.</param>
        /// <param name="otherCoefficient">A scalar that multiplies each entry of <paramref name="otherVector"/>.</param>
        public void AxpyIntoThis(Vector2 otherVector, double otherCoefficient)
        {
            this.data[0] += otherCoefficient * otherVector.data[0];
            this.data[1] += otherCoefficient * otherVector.data[1];
        }

        /// <summary>
        /// See <see cref="IVector.AxpySubvectorIntoThis(int, IVectorView, double, int, int)"/>
        /// </summary>
        public void AxpySubvectorIntoThis(int destinationIndex, IVectorView sourceVector, double sourceCoefficient, 
            int sourceIndex, int length)
        {
            Preconditions.CheckSubvectorDimensions(this, destinationIndex, length);
            Preconditions.CheckSubvectorDimensions(sourceVector, sourceIndex, length);

            if (sourceVector is Vector2 casted)
            {
                for (int i = 0; i < length; ++i) data[i + destinationIndex] += sourceCoefficient * casted.data[i + sourceIndex];
            }
            else
            {
                for (int i = 0; i < length; ++i) data[i + destinationIndex] += sourceCoefficient * sourceVector[i + sourceIndex];
            }
        }

        /// <summary>
        /// See <see cref="IVector.Clear"/>.
        /// </summary>
        public void Clear() //TODO: Is Array.Clear faster here?
        {
            data[0] = 0.0;
            data[1] = 0.0;
        }

        /// <summary>
        /// See <see cref="IVector.Copy(bool)"/>
        /// </summary>
        public IVector Copy(bool copyIndexingData = false)
            => Copy();

        /// <summary>
        /// Initializes a new instance of <see cref="Vector2"/> by deep copying the entries of this instance.
        /// </summary>
        public Vector2 Copy() => new Vector2(new double[] { data[0], data[1] });
        

        /// <summary>
        /// See <see cref="IVector.CopyFrom(IVectorView)"/>
        /// </summary>
        public void CopyFrom(IVectorView sourceVector)
        {
            Preconditions.CheckVectorDimensions(this, sourceVector);
            if (sourceVector is Vector2 casted)
            {
                data[0] = casted.data[0];
                data[1] = casted.data[1];
            }
            else
            {
                data[0] = sourceVector[0];
                data[1] = sourceVector[1];
            }
        }

        /// <summary>
        /// See <see cref="IVector.CopyNonContiguouslyFrom(IVectorView, int[])"/>
        /// </summary>
        public void CopyNonContiguouslyFrom(IVectorView otherVector, int[] otherIndices)
            => DenseStrategies.CopyNonContiguouslyFrom(this, otherVector, otherIndices);

        /// <summary>
        /// See <see cref="IVector.CopySubvectorFrom(int, IVectorView, int, int)"/>
        /// </summary>
        public void CopySubvectorFrom(int destinationIndex, IVectorView sourceVector, int sourceIndex, int length)
        {
            Preconditions.CheckSubvectorDimensions(this, destinationIndex, length);
            Preconditions.CheckSubvectorDimensions(sourceVector, sourceIndex, length);

            if (sourceVector is Vector2 casted)
            {
                for (int i = 0; i < length; ++i) data[i + destinationIndex] = casted.data[i + sourceIndex];
            }
            else
            {
                for (int i = 0; i < length; ++i) data[i + destinationIndex] = sourceVector[i + sourceIndex];
            }
        }

        /// <summary>
        /// See <see cref="IVector.CopyNonContiguouslyFrom(int[], IVectorView, int[])"/>
        /// </summary>
        public void CopyNonContiguouslyFrom(int[] thisIndices, IVectorView otherVector, int[] otherIndices)
            => DenseStrategies.CopyNonContiguouslyFrom(this, thisIndices, otherVector, otherIndices);

        /// <summary>
        /// See <see cref="IVectorView.CopyToArray"/>.
        /// </summary>
        public double[] CopyToArray() => new double[] { data[0], data[1] };

        /// <summary>
        /// See <see cref="IVectorView.CreateZeroVectorWithSameFormat"/>
        /// </summary>
        public IVector CreateZeroVectorWithSameFormat() => new Vector2(new double[2]);

        /// <summary>
        /// Performs the operation: result = this[0] * other[1] - this[1] * other[0]. The result is a scalar value.  
        /// Also note that: other.Cross(this) = - this.Cross(other).
        /// </summary>
        /// <param name="vector">A vector with two entries.</param>
        public double CrossProduct(Vector2 vector) => this.data[0] * vector.data[1] - this.data[1] * vector.data[0];

        /// <summary>
        /// See <see cref="IEntrywiseOperableView1D{TVectorIn, TVectorOut}.DoEntrywise(TVectorIn, Func{double, double, double})"/>.
        /// </summary>
        public IVector DoEntrywise(IVectorView vector, Func<double, double, double> binaryOperation)
        {
            if (vector is Vector2 casted) return DoEntrywise(vector, binaryOperation);
            else
            {
                Preconditions.CheckVectorDimensions(this, vector);
                return new Vector2(new double[] 
                {
                    binaryOperation(this.data[0], vector[0]),
                    binaryOperation(this.data[1], vector[1])
                });
            }
        }

        /// <summary>
        /// See <see cref="IEntrywiseOperableView1D{TVectorIn, TVectorOut}.DoEntrywise(TVectorIn, Func{double, double, double})"/>.
        /// </summary>
        public Vector2 DoEntrywise(Vector2 vector, Func<double, double, double> binaryOperation)
        {
            return new Vector2(new double[] { binaryOperation(this.data[0], vector.data[0]),
                binaryOperation(this.data[1], vector.data[1]) });
        }

        /// <summary>
        /// See <see cref="IEntrywiseOperable1D{TVectorIn}.DoEntrywiseIntoThis(TVectorIn, Func{double, double, double})"/>
        /// </summary>
        public void DoEntrywiseIntoThis(IVectorView otherVector, Func<double, double, double> binaryOperation)
        {
            if (otherVector is Vector2 casted) DoEntrywiseIntoThis(casted, binaryOperation);
            else
            {
                Preconditions.CheckVectorDimensions(this, otherVector);
                data[0] = binaryOperation(data[0], otherVector[0]);
                data[1] = binaryOperation(data[1], otherVector[1]);
            }
        }

        /// <summary>
        /// See <see cref="IEntrywiseOperable1D{TVectorIn}.DoEntrywiseIntoThis(TVectorIn, Func{double, double, double})"/>
        /// </summary>
        public void DoEntrywiseIntoThis(Vector2 vector, Func<double, double, double> binaryOperation)
        {
            this.data[0] = binaryOperation(this.data[0], vector.data[0]);
            this.data[1] = binaryOperation(this.data[1], vector.data[1]);
        }

        /// <summary>
        /// See <see cref="IEntrywiseOperableView1D{TVectorIn, TVectorOut}.DoToAllEntries(Func{double, double})"/>.
        /// </summary>
        IVector IEntrywiseOperableView1D<IVectorView, IVector>.DoToAllEntries(Func<double, double> unaryOperation) 
            => DoToAllEntries(unaryOperation);

        /// <summary>
        /// See <see cref="IEntrywiseOperableView1D{TVectorIn, TVectorOut}.DoToAllEntries(Func{double, double})"/>.
        /// </summary>
        public Vector2 DoToAllEntries(Func<double, double> unaryOperation)
        {
            return new Vector2(new double[] { unaryOperation(data[0]), unaryOperation(data[1]) });
        }

        /// <summary>
        /// See <see cref="IEntrywiseOperable1D{TVectorIn}.DoToAllEntriesIntoThis(Func{double, double})"/>
        /// </summary>
        public void DoToAllEntriesIntoThis(Func<double, double> unaryOperation)
        {
            this.data[0] = unaryOperation(this.data[0]);
            this.data[1] = unaryOperation(this.data[1]);
        }

        /// <summary>
        /// See <see cref="IVectorView.DotProduct(IVectorView)"/>.
        /// </summary>
        public double DotProduct(IVectorView vector)
        {
            if (vector is Vector2 casted) return DotProduct(casted);
            else
            {
                Preconditions.CheckVectorDimensions(this, vector);
                return data[0] * vector[0] + data[1] * vector[1];
            }
        }

        /// <summary>
        /// Calculates the dot (or inner/scalar) product of this vector with <paramref name="vector"/>:
        /// result = this[0] * <paramref name="vector"/>[0] + this[1] * <paramref name="vector"/>[1].
        /// </summary>
        /// <param name="vector">A vector with two entries</param>
        public double DotProduct(Vector2 vector) => this.data[0] * vector.data[0] + this.data[1] * vector.data[1];

        /// <summary>
        /// See <see cref="IIndexable1D.Equals(IIndexable1D, double)"/>.
        /// </summary>
        bool IIndexable1D.Equals(IIndexable1D other, double tolerance)
        {
            if (other.Length != 2) return false;
            else
            {
                var comparer = new ValueComparer(tolerance);
                return comparer.AreEqual(this.data[0], other[0]) && comparer.AreEqual(this.data[1], other[1]);
            }
        }

        /// <summary>
        /// Returns true if this[i] - <paramref name="other"/>[i] is within the acceptable <paramref name="tolerance"/> for all
        /// 0 &lt;= i &lt; 2. 
        /// </summary>
        /// <param name="other">The other vector that this <see cref="Vector2"/> instance will be compared to.</param>
        /// <param name="tolerance">The entries at index i of the two vectors will be considered equal, if
        ///     (<paramref name="other"/>[i] - this[i]) / this[i] &lt;= <paramref name="tolerance"/>. Setting 
        ///     <paramref name="tolerance"/> = 0, will check if these entries are exactly the same.</param>
        public bool Equals(Vector2 other, double tolerance = 1e-13)
        {
            var comparer = new ValueComparer(tolerance);
            return comparer.AreEqual(this.data[0], other.data[0]) && comparer.AreEqual(this.data[1], other.data[1]);
        }

        /// <summary>
        /// See <see cref="IVectorView.LinearCombination(double, IVectorView, double)"/>.
        /// </summary>
        public IVector LinearCombination(double thisCoefficient, IVectorView otherVector, double otherCoefficient)
        {
            if (otherVector is Vector2 casted) return LinearCombination(thisCoefficient, casted, otherCoefficient);
            else if (thisCoefficient == 1.0) return Axpy(otherVector, otherCoefficient);
            else
            {
                Preconditions.CheckVectorDimensions(this, otherVector);
                return new Vector2(new double[]
                {
                    thisCoefficient * data[0] + otherCoefficient * otherVector[0],
                    thisCoefficient * data[1] + otherCoefficient * otherVector[1]
                });
            }
        }

        /// <summary>
        /// Performs the operation: result[i] = <paramref name="thisCoefficient"/> * this[i] + 
        ///     <paramref name="otherCoefficient"/> * <paramref name="otherVector"/>[i], 
        /// for 0 &lt;= i &lt; 2. The resulting vector is written to a new <see cref="Vector2"/> and then returned.
        /// </summary>
        /// <param name="thisCoefficient">A scalar that multiplies each entry of this vector.</param>
        /// <param name="otherVector">A vector with two entries.</param>
        /// <param name="otherCoefficient">A scalar that multiplies each entry of <paramref name="otherVector"/>.</param>
        public Vector2 LinearCombination(double thisCoefficient, Vector2 otherVector, double otherCoefficient)
        {
            if (thisCoefficient == 1.0) return Axpy(otherVector, otherCoefficient);
            else return new Vector2(new double[] 
            {
                thisCoefficient * this.data[0] + otherCoefficient * otherVector.data[0],
                thisCoefficient * this.data[1] + otherCoefficient * otherVector.data[1]
            });
        }

        /// <summary>
        /// See <see cref="IVector.LinearCombinationIntoThis(double, IVectorView, double)"/>
        /// </summary>
        public void LinearCombinationIntoThis(double thisCoefficient, IVectorView otherVector, double otherCoefficient)
        {
            if (otherVector is Vector2 casted) LinearCombinationIntoThis(thisCoefficient, casted, otherCoefficient);
            else if (thisCoefficient == 1.0) AxpyIntoThis(otherVector, otherCoefficient);
            else
            {
                
                Preconditions.CheckVectorDimensions(this, otherVector);
                data[0] = thisCoefficient * data[0] * otherCoefficient * otherVector[0];
                data[1] = thisCoefficient * data[1] * otherCoefficient * otherVector[1];
            }
        }

        /// <summary>
        /// Performs the operation: this[i] = <paramref name="thisCoefficient"/> * this[i] + 
        ///     <paramref name="otherCoefficient"/> * <paramref name="otherVector"/>[i], 
        /// for 0 &lt;= i &lt; 2. The resulting vector overwrites the entries of this <see cref="Vector2"/> instance.
        /// </summary>
        /// <param name="thisCoefficient">A scalar that multiplies each entry of this vector.</param>
        /// <param name="otherVector">A vector with two entries.</param>
        /// <param name="otherCoefficient">A scalar that multiplies each entry of <paramref name="otherVector"/>.</param>
        public void LinearCombinationIntoThis(double thisCoefficient, Vector2 otherVector, double otherCoefficient)
        {
            if (thisCoefficient == 1.0) AxpyIntoThis(otherVector, otherCoefficient);
            else
            {
                this.data[0] = thisCoefficient * this.data[0] * otherCoefficient * otherVector.data[0];
                this.data[1] = thisCoefficient * this.data[1] * otherCoefficient * otherVector.data[1];
            }
        }

        /// <summary>
        /// See <see cref="IVectorView.Norm2"/>
        /// </summary>
        public double Norm2() => Math.Sqrt(data[0] * data[0] + data[1] * data[1]);

        /// <summary>
        /// See <see cref="IReducible.Reduce(double, ProcessEntry, ProcessZeros, Reduction.Finalize)"/>.
        /// </summary>
        public double Reduce(double identityValue, ProcessEntry processEntry, ProcessZeros processZeros, Finalize finalize)
        {
            double accumulator = identityValue; // no zeros implied
            accumulator = processEntry(data[0], accumulator);
            accumulator = processEntry(data[1], accumulator);
            return finalize(accumulator);
        }

        /// <summary>
        /// See <see cref="IVectorView.Scale(double)"/>.
        /// </summary>
        IVector IVectorView.Scale(double scalar) => Scale(scalar);

        /// <summary>
        /// Performs the operation: result[i] = <paramref name="scalar"/> * this[i],
        /// for 0 &lt;= i &lt; 2. The resulting vector is written to a new <see cref="Vector2"/> and then returned.
        /// </summary>
        /// <param name="scalar">A scalar that multiplies each entry of this vector.</param>
        public Vector2 Scale(double scalar) => new Vector2(new double[] { scalar * data[0], scalar * data[1] });

        /// <summary>
        /// Performs the operation: this[i] = <paramref name="scalar"/> * this[i],
        /// for 0 &lt;= i &lt; 2. The resulting vector overwrites the entries of this <see cref="Vector2"/> instance.
        /// </summary>
        /// <param name="scalar">A scalar that multiplies each entry of this vector.</param>
        public void ScaleIntoThis(double scalar)
        {
            data[0] *= scalar;
            data[1] *= scalar;
        }

        /// <summary>
        /// See <see cref="IVector.Set(int, double)"/>
        /// </summary>
        public void Set(int index, double value)
        {
            data[index] = value;
        }

        /// <summary>
        /// Sets all entries of this vector to be equal to <paramref name="value"/>.
        /// </summary>
        /// <param name="value">The value that all entries of the this vector will be equal to.</param>
        public void SetAll(double value)
        {
            data[0] = value;
            data[1] = value;
        }

        /// <summary>
        /// Performs the operation: this[i] = this[i] - <paramref name="vector"/>[i], 
        /// for 0 &lt;= i &lt; 2. The resulting vector overwrites the entries of this <see cref="Vector2"/> instance.
        /// </summary>
        /// <param name="vector">A vector with two entries.</param>
        public void SubtractIntoThis(Vector2 vector)
        {
            this.data[0] -= vector.data[0];
            this.data[1] -= vector.data[1];
        }
    }
}
