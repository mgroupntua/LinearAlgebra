//TODO: Use a dedicated DiagonalMatrix class, instead of passing in double[] or Vector. It will implement the inverse and 
//      multiplication routines. It will also handle distributed matrices. E.g. IDiagonal IMatrixView.GetDiagonal() which will 
//      then have an IDiagonalMatrix.Inverse(). The problem is how we will go from CSR to DiagonalMatrix. Perhaps it would be 
//      better to use the DOK instead.
//TODO: Alternative: instead of demanding the caller to extract the diagonal, this class should read the matrix and only access 
//      its diagonal. I think this alternative is less flexible and more difficult to implement.
namespace MGroup.LinearAlgebra.Iterative.Preconditioning
{
	using System;

	using MGroup.LinearAlgebra.Commons;
	using MGroup.LinearAlgebra.Exceptions;
	using MGroup.LinearAlgebra.Matrices;
	using MGroup.LinearAlgebra.Vectors;

	/// <summary>
	/// Implements the Jacobi or diagonal preconditioner for a square matrix. If A is the original matrix, the Jacobi  
	/// preconditioner is a matrix M, such that it oncly contains the diagonal of A and inverse(M) is also diagonal with 
	/// entries: 1/A[0,0], 1/A[1,1], ... The Jacobi preconditioner is cheapest to build and apply, but doesn't improve 
	/// convergence as much as other preconditioners.
	/// </summary>
	public class JacobiPreconditioner : IPreconditioner
	{
		public const double DefaultTolerance = 1e-10;
		private readonly double tolerance;
		private double[] inverseDiagonal;

		/// <summary>
		/// Initializes a new instance of <see cref="JacobiPreconditioner"/> for the linear system's matrix whose main diagonal
		/// is provided in <paramref name="diagonal"/>.
		/// </summary>
		/// <param name="tolerance">The value under which a diagonal entry will be considered as zero.</param>
		/// <exception cref="SingularMatrixException">If there is a zero diagonal entry.</exception>
		public JacobiPreconditioner(double tolerance = DefaultTolerance)
		{
			this.tolerance = tolerance;
		}

		public IPreconditioner CopyWithInitialSettings() => new JacobiPreconditioner(tolerance);

		/// <summary>
		/// <inheritdoc/>
		/// </summary>
		public void SolveLinearSystem(IReadOnlyVector rhsVector, IVector lhsVector)
		{
			Preconditions.CheckSystemSolutionDimensions(inverseDiagonal.Length, rhsVector.Length);
			if ((rhsVector is Vector rhs) && (lhsVector is Vector lhs))
			{
				double[] y = rhs.RawData;
				double[] x = lhs.RawData;
				for (int i = 0; i < inverseDiagonal.Length; ++i)
				{
					x[i] = inverseDiagonal[i] * y[i];
				}
			}
			else
			{
				for (int i = 0; i < inverseDiagonal.Length; ++i)
				{
					lhsVector.Set(i, inverseDiagonal[i] * rhsVector[i]);
				}
			}
		}

		public void UpdateMatrix(IReadOnlyMatrix matrix, bool isPatternModified)
		{
			this.inverseDiagonal = matrix.GetDiagonalAsArray();
			for (int i = 0; i < inverseDiagonal.Length; ++i)
			{
				double val = inverseDiagonal[i];
				if (Math.Abs(val) <= tolerance)
				{
					throw new SingularMatrixException($"Zero diagonal entry at index {i}");
				}

				this.inverseDiagonal[i] = 1.0 / val;
			}
		}
	}
}
