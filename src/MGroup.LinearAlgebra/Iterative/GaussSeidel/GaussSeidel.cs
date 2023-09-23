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
		private readonly ISolutionConvergenceCriterion convergenceCriterion;
		private readonly double convergenceTolerance;
		private readonly bool forwardGaussSeidel;
		private readonly IMaxIterationsProvider maxIterationsProvider;

		public GaussSeidel(ISolutionConvergenceCriterion convergenceCriterion, double convergenceTolerance, bool forwardGaussSeidel, IMaxIterationsProvider maxIterationsProvider)
		{
			this.convergenceCriterion = convergenceCriterion;
			this.convergenceTolerance = convergenceTolerance;
			this.forwardGaussSeidel = forwardGaussSeidel;
			this.maxIterationsProvider = maxIterationsProvider;
		}

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
		/// <param name="initialGuessIsZero">
		/// If <paramref name="solution"/> is 0, then set <paramref name="initialGuessIsZero"/> to true to avoid performing the
		/// operation b-A*0 before starting.
		/// </param>
		/// <exception cref="NonMatchingDimensionsException">
		/// Thrown if <paramref name="rhs"/> or <paramref name="solution"/> violate the described constraints.
		/// </exception>
		public IterativeStatistics Solve(CsrMatrix matrix, Vector rhs, Vector solution, bool initialGuessIsZero = false)
		{
			Preconditions.CheckSquareLinearSystemDimensions(matrix, solution, rhs);

			int n = matrix.NumColumns;
			int maxIterations = maxIterationsProvider.GetMaxIterations(n);
			var previousSolution = Vector.CreateZero(n);
			double convergenceMetric = double.MaxValue;
			int iter = 0;
			while (iter < maxIterations)
			{
				previousSolution.CopyFrom(solution);
				if (forwardGaussSeidel)
				{
					CsrMultiplications.GaussSeidelForwardIteration(n, matrix.RawValues, matrix.RawRowOffsets, matrix.RawColIndices, solution.RawData, rhs.RawData);
				}
				else
				{
					CsrMultiplications.GaussSeidelBackwardIteration(n, matrix.RawValues, matrix.RawRowOffsets, matrix.RawColIndices, solution.RawData, rhs.RawData);
				}
				++iter; // Each algorithm iteration corresponds to one matrix-vector multiplication or, in this case, GS iteration

				convergenceMetric = convergenceCriterion.CalculateConvergenceMetric(solution, previousSolution);
				if (convergenceMetric < convergenceTolerance)
				{
					break;
				}
			}

			return new IterativeStatistics
			{
				AlgorithmName = name,
				HasConverged = iter < maxIterations,
				NumIterationsRequired = iter,
				ConvergenceCriterion = (convergenceCriterion.DescribeConvergenceCriterion(convergenceTolerance), convergenceMetric),
			};
		}

		/// <summary>
		/// Constructs <see cref="GaussSeidel"/> instances, allows the user to specify some or all of the required parameters and 
		/// provides defaults for the rest.
		/// Author: Constantinos Atzarakis
		/// </summary>
		public class Builder
		{
			public ISolutionConvergenceCriterion ConvergenceCriterion { get; set; } = new AbsoluteSolutionConvergenceCriterion();

			public double ConvergenceTolerance { get; set; } = 1E-10;


			public bool ForwardGaussSeidel { get; set; } = true;

			/// <summary>
			/// Specifies how to calculate the maximum iterations that the GS algorithm will run for.
			/// </summary>
			public IMaxIterationsProvider MaxIterationsProvider { get; set; } = new PercentageMaxIterationsProvider(1.0);


			/// <summary>
			/// Creates a new instance of <see cref="GaussSeidel"/>.
			/// </summary>
			public GaussSeidel Build()
				=> new GaussSeidel(ConvergenceCriterion, ConvergenceTolerance, ForwardGaussSeidel, MaxIterationsProvider);
		}
	}
}
