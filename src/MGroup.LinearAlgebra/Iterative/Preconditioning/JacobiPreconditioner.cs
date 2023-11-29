using System;
using MGroup.LinearAlgebra.Commons;
using MGroup.LinearAlgebra.Exceptions;
using MGroup.LinearAlgebra.Matrices;
using MGroup.LinearAlgebra.Vectors;

//TODO: Use a dedicated DiagonalMatrix class, instead of passing in double[] or Vector. It will implement the inverse and 
//      multiplication routines. It will also handle distributed matrices. E.g. IDiagonal IMatrixView.GetDiagonal() which will 
//      then have an IDiagonalMatrix.Inverse(). The problem is how we will go from CSR to DiagonalMatrix. Perhaps it would be 
//      better to use the DOK instead.
//TODO: Alternative: instead of demanding the caller to extract the diagonal, this class should read the matrix and only access 
//      its diagonal. I think this alternative is less flexible and more difficult to implement.
namespace MGroup.LinearAlgebra.Iterative.Preconditioning
{
    /// <summary>
    /// Implements the Jacobi or diagonal preconditioner for a square matrix. If A is the original matrix, the Jacobi  
    /// preconditioner is a matrix M, such that it oncly contains the diagonal of A and inverse(M) is also diagonal with 
    /// entries: 1/A[0,0], 1/A[1,1], ... The Jacobi preconditioner is cheapest to build and apply, but doesn't improve 
    /// convergence as much as other preconditioners.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class JacobiPreconditioner: IPreconditioner
    {
        public const double DefaultTolerance = 1e-10;
		private readonly bool isDiagonalInverted;
        private readonly double[] diagonal;

		/// <summary>
		/// Initializes a new instance of <see cref="JacobiPreconditioner"/> for the linear system's matrix whose main diagonal
		/// is provided in <paramref name="diagonal"/>.
		/// </summary>
		/// <param name="diagonal">
		/// The main diagonal of the original matrix of the linear system. Constraints: all its entries must be non-zero.
		/// </param>
		/// <param name="preinvert">
		/// If true, the diagonal will be be iverted at this point and <see cref="SolveLinearSystem(IVectorView, IVector)"/> will
		/// multiply with the values of the inverse diagonal, which is faster overall. If false, 
		/// <see cref="SolveLinearSystem(IVectorView, IVector)"/> will divide with the values of the original diagonal, which is
		/// slower, but may be more stable numerically, when the values are very small.
		/// </param>
		/// <param name="tolerance">
		/// The value under which a diagonal entry will be considered as zero. Will not be used if 
		/// <paramref name="preinvert"/> == false.
		/// </param>
		/// <exception cref="SingularMatrixException">If there is a zero diagonal entry.</exception>
		public JacobiPreconditioner(double[] diagonal, bool preinvert = true, double tolerance = DefaultTolerance)
        {
            Order = diagonal.Length;
            this.diagonal = new double[Order];
			isDiagonalInverted = preinvert;
			if (preinvert)
			{
				for (int i = 0; i < Order; ++i)
				{
					double val = diagonal[i];
					if (Math.Abs(val) <= tolerance) throw new SingularMatrixException($"Zero diagonal entry at index {i}");
					this.diagonal[i] = 1.0 / val;
				}
			}
			else
			{
				this.diagonal = diagonal;
			}
        }

        /// <summary>
        /// The number of rows/columns of the preconditioner and the original matrix
        /// </summary>
        public int Order { get; }

        /// <summary>
        /// See <see cref="IPreconditioner.SolveLinearSystem(Vector)"/>
        /// </summary>
        public void SolveLinearSystem(IVectorView rhsVector, IVector lhsVector)
        {
            Preconditions.CheckSystemSolutionDimensions(Order, rhsVector.Length);
			if (isDiagonalInverted)
			{
				for (int i = 0; i < Order; ++i)
				{
					lhsVector.Set(i, diagonal[i] * rhsVector[i]);
				}
			}
            else
			{
				for (int i = 0; i < Order; ++i)
				{
					lhsVector.Set(i, rhsVector[i] / diagonal[i]);
				}
			}
        }

        /// <summary>
        /// Creates instances of <see cref="JacobiPreconditioner"/>.
        /// </summary>
        public class Factory: IPreconditionerFactory
        {
            private readonly double tolerance;

			/// <summary>
			/// Initializes a new instance of <see cref="JacobiPreconditioner.Factory"/> with the specified settings.
			/// </summary>
			/// <param name="tolerance">The value under which a diagonal entry will be considered as zero.</param>
			public Factory() { }

			/// <summary>
			/// The value under which a diagonal entry will be considered as zero. Will not be used if 
			/// <see cref="PreInvert"/> == false.
			/// </summary>
			public double InversionTolerance { get; set; } = JacobiPreconditioner.DefaultTolerance;

			/// <summary>
			/// If true, the diagonal will be be iverted at this point and <see cref="SolveLinearSystem(IVectorView, IVector)"/> 
			/// will multiply with the values of the inverse diagonal, which is faster overall. If false, 
			/// <see cref="SolveLinearSystem(IVectorView, IVector)"/> will divide with the values of the original diagonal, 
			/// which is slower, but may be more stable numerically, when the values are very small.
			/// </summary>
			public bool PreInvert { get; set; } = true;

            /// <summary>
            /// See <see cref="IPreconditionerFactory.CreatePreconditionerFor(IMatrixView)"/>.
            /// </summary>
            public IPreconditioner CreatePreconditionerFor(IMatrixView matrix) 
                => new JacobiPreconditioner(matrix.GetDiagonalAsArray(), PreInvert, InversionTolerance);
        }
    }
}
