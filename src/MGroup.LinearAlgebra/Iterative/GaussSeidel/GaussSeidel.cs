using System;

using MGroup.LinearAlgebra.Exceptions;
using MGroup.LinearAlgebra.Iterative.Termination;
using MGroup.LinearAlgebra.Matrices;
using MGroup.LinearAlgebra.Vectors;
using MGroup.LinearAlgebra.Commons;

namespace MGroup.LinearAlgebra.Iterative
{
	using Reduction;

	/// <summary>
    /// Implements the Gauss-Seidel algorithm for solving linear systems.
    /// Convergence is guaranteed only for strictly diagonally dominant or positive definite (symmetric) matrices.
    /// Might converge in general matrix systems, as well, but with no guarantee.
    /// Authors: Constantinos Atzarakis
    /// </summary>
	public class GaussSeidel
    {
        private const string name = "Gauss-Seidel";
        private readonly IMaxIterationsProvider maxIterationsProvider;
        private readonly double residualTolerance;

        private IVector residual;
        private IVector solution;

		public GaussSeidel(double residualTolerance, IMaxIterationsProvider maxIterationsProvider)
        {
			this.maxIterationsProvider = maxIterationsProvider;
			this.residualTolerance = residualTolerance;
		}


        /// <summary>
        /// The current iteration of the algorithm. It belongs to the interval [0, maxIterations).
        /// </summary>
        public int Iteration { get; private set; }

        /// <summary>
        /// The matrix A of the linear system or another object that implements matrix-vector multiplications.
        /// </summary>
        public ILinearTransformation Matrix { get; private set; }

		/// <summary>
        /// The residual vector r = b - A * x.
        /// </summary>
        public IVectorView Residual => residual;

        /// <summary>
        /// The right hand side of the linear system b = A * x.
        /// </summary>
        public IVectorView Rhs { get; private set; }

        /// <summary>
        /// The current approximation to the solution of the linear system A * x = b
        /// </summary>
        public IVectorView Solution => solution;

        /// <summary>
        /// Releases references to the vectors and matrices used by this object and sets scalars to their default values.
        /// </summary>
        public void Clear()
        {
	        // throw new NotImplementedException();
	        Matrix = null;
	        Rhs = null;
	        solution = null;
	        residual = null;
	        Iteration = -1;
        }

        public IterativeStatistics Solve(IMatrixView matrix, IVectorView rhs, IVector solution)
            => Solve(new ExplicitMatrixTransformation(matrix), rhs, solution);

        /// <summary>
        /// Solves the linear system A * x = b, where A = <paramref name="matrix"/> and b = <paramref name="rhs"/>.
        /// Initially x = <paramref name="initialGuess"/> and then it converges to the solution.
        /// </summary>
        /// <param name="matrix">The matrix A of the linear system A * x = b. It must be symmetric positive definite
        /// or strictly diagonally dominant for ensured convergence.</param>
        /// <param name="rhs">
        /// The right hand side vector b of the linear system A * x = b. Constraints:
        /// <paramref name="rhs"/>.<see cref="IIndexable1D.Length"/> 
        /// == <paramref name="matrix"/>.<see cref="IIndexable2D.NumRows"/>.
        /// </param>
        /// <param name="solution">
        /// The vector from which to start refining the solution vector x. Constraints:
        /// <paramref name="solution"/>.<see cref="IIndexable1D.Length"/>
        /// == <paramref name="matrix"/>.<see cref="IIndexable2D.NumColumns"/>.
        /// </param>
        /// <exception cref="NonMatchingDimensionsException">
        /// Thrown if <paramref name="rhs"/> or <paramref name="solution"/> violate the described constraints.
        /// </exception>
        public IterativeStatistics Solve(ILinearTransformation matrix, IVectorView rhs, IVector solution)
        {
            // throw new NotImplementedException();
             //TODO: these will also be checked by the matrix vector multiplication.
             Preconditions.CheckMultiplicationDimensions(matrix.NumColumns, solution.Length);
             Preconditions.CheckSystemSolutionDimensions(matrix.NumRows, rhs.Length);

             this.Matrix = matrix;
             this.Rhs = rhs;
             this.solution = solution;

             // r = b - A * x
             residual = ExactResidual.Calculate(matrix, rhs, solution);

             return SolveInternal(maxIterationsProvider.GetMaxIterations(matrix.NumColumns));
        }

        /// <summary>
        /// Main body of Gauss-Seidel algorithm <see cref="GaussSeidel"/>
        /// Author: Constantinos Atzarakis
        /// </summary>
        private IterativeStatistics SolveInternal(int maxIterations)
        {
			throw new NotImplementedException();
		}

        /// <summary>
        /// Constructs <see cref="GaussSeidel"/> instances, allows the user to specify some or all of the required parameters and 
        /// provides defaults for the rest.
        /// Author: Constantinos Atzarakis
        /// </summary>
        public class Builder
        {
            /// <summary>
            /// Specifies how to calculate the maximum iterations that the GS algorithm will run for.
            /// </summary>
            public IMaxIterationsProvider MaxIterationsProvider { get; set; } = new PercentageMaxIterationsProvider(1.0);

            public double ResidualTolerance { get; set; } = 1E-10;

			/// <summary>
			/// Creates a new instance of <see cref="GaussSeidel"/>.
			/// </summary>
            public GaussSeidel Build()
                => new GaussSeidel(ResidualTolerance, MaxIterationsProvider);
        }
    }
}
