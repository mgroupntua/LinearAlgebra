using System;
using System.Collections.Generic;
using MGroup.LinearAlgebra.Commons;
using MGroup.LinearAlgebra.Exceptions;
using MGroup.LinearAlgebra.Matrices;
using MGroup.LinearAlgebra.Output.Formatting;
using MGroup.LinearAlgebra.Providers;
using MGroup.LinearAlgebra.Vectors;

namespace MGroup.LinearAlgebra.Triangulation
{
    /// <summary>
    /// Base class for factorizations (triangulations in particular) where the matrix is stored in Skyline format.
    /// </summary>
    public abstract class SkylineFactorizationBase: IIndexable2D, ISparseMatrix, ITriangulation
    {
        protected readonly double[] values;
        protected readonly int[] diagOffsets;

        protected SkylineFactorizationBase(int order, double[] values, int[] diagOffsets)
        {
            this.Order = order;
            this.values = values;
            this.diagOffsets = diagOffsets;
        }

		/// <summary>
		/// Matrix symmetry properties (by default symmetric).
		/// </summary>
		public MatrixSymmetry MatrixSymmetry => MatrixSymmetry.Symmetric;

		/// <summary>
		/// The number of columns of the matrix. 
		/// </summary>
		public int NumColumns => Order;

        /// <summary>
        /// The number of rows of the matrix.
        /// </summary>
        public int NumRows => Order;

        public int Order { get; }

		/// <summary>
		/// The internal array that stores the non-zero entries of the matrix's upper triangle in column major order, 
		/// starting from the diagonal and going upwards. Its length is equal to the number of non-zero entries. 
		/// </summary>
		public double[] RawValues => values;

		/// <summary>
		/// The internal array that stores the indices into <see cref="RawValues"/> of the diagonal entries of the matrix. 
		/// Its length = order + 1, with the last entry being equal to nnz.
		/// </summary>
		public int[] RawDiagOffsets => diagOffsets;

		/// <summary>
		/// See <see cref="IIndexable2D.this[int, int]"/>.
		/// </summary>
		public double this[int rowIdx, int colIdx]
        {
            get
            {
                SkylineMatrix.ProcessIndices(ref rowIdx, ref colIdx);
                int entryHeight = colIdx - rowIdx; // excluding diagonal
                int maxColumnHeight = diagOffsets[colIdx + 1] - diagOffsets[colIdx] - 1; // excluding diagonal
                if (entryHeight > maxColumnHeight) return 0.0; // outside stored non zero pattern
                else return values[diagOffsets[colIdx] + entryHeight];
            }
        }

        /// <summary>
        /// See <see cref="ISparseMatrix.CountNonZeros"/>.
        /// </summary>
        public int CountNonZeros() => values.Length;

        /// <summary>
        /// See <see cref="ISparseMatrix.EnumerateNonZeros"/>.
        /// </summary>
        public IEnumerable<(int row, int col, double value)> EnumerateNonZeros()
            => SkylineMatrix.CreateFromArrays(Order, values, diagOffsets, false, false).EnumerateNonZeros();

        /// <summary>
        /// See <see cref="IIndexable2D.Equals(IIndexable2D, double)"/>.
        /// </summary>
        public bool Equals(IIndexable2D other, double tolerance = 1E-13) //TODO: what are the semantics of this? It cannot be compared to matrices. Perhaps IIndexable2D should not have Equals()
            => SkylineMatrix.CreateFromArrays(Order, values, diagOffsets, false, false).Equals(other, tolerance);

        /// <summary>
        /// See <see cref="ISparseMatrix.GetSparseFormat"/>.
        /// </summary>
        public SparseFormat GetSparseFormat()
            => SkylineMatrix.CreateFromArrays(Order, values, diagOffsets, false, false).GetSparseFormat();


        public abstract double CalcDeterminant();
        public abstract void SolveLinearSystem(Vector rhs, Vector solution);
    }
}
