namespace MGroup.LinearAlgebra.Iterative.Stationary
{
	using MGroup.LinearAlgebra.Commons;
	using MGroup.LinearAlgebra.Exceptions;
	using MGroup.LinearAlgebra.Iterative.Termination.Convergence;
	using MGroup.LinearAlgebra.Iterative.Termination.Iterations;
	using MGroup.LinearAlgebra.Matrices;
	using MGroup.LinearAlgebra.Vectors;

	/// <summary>
	/// Implements a general stationary iterative method algorithm for solving linear systems.
	/// Convergence depends on the specific algorithm used.
	/// </summary>
	public class StationaryIterativeMethod
	{
		private readonly ISolutionConvergenceCriterion convergenceCriterion;
		private readonly double convergenceTolerance;
		private readonly IStationaryIteration stationaryIteration;
		private readonly IMaxIterationsProvider maxIterationsProvider;

		public StationaryIterativeMethod(IStationaryIteration stationaryIteration, double convergenceTolerance, 
			ISolutionConvergenceCriterion convergenceCriterion, IMaxIterationsProvider maxIterationsProvider)
		{
			this.stationaryIteration = stationaryIteration;
			this.convergenceCriterion = convergenceCriterion;
			this.convergenceTolerance = convergenceTolerance;
			this.maxIterationsProvider = maxIterationsProvider;
		}

		/// <summary>
		/// Solves the linear system A * x = b, where A = <paramref name="matrix"/> and b = <paramref name="rhs"/>.
		/// Initially x = <paramref name="initialGuess"/> and then it converges to the solution.
		/// </summary>
		/// <param name="matrix">
		/// The matrix A of the linear system A * x = b. A must be symmetric positive definite or strictly diagonally dominant 
		/// for ensured convergence.
		/// </param>
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
		public IterativeStatistics Solve(IReadOnlyMatrix matrix, IReadOnlyVector rhs, IVector solution)
		{
			Preconditions.CheckSquareLinearSystemDimensions(matrix, solution, rhs);

			stationaryIteration.UpdateMatrix(matrix, true);
			Vector rhsDense = (Vector)rhs;
			Vector solutionDense = (Vector)solution;
			Vector previousSolution = solutionDense.Copy();

			double convergenceMetric = double.NaN;
			int maxIterations = maxIterationsProvider.GetMaxIterations(matrix.NumRows);
			int iter = 0;

			while (iter < maxIterations)
			{
				previousSolution.CopyFrom(solution);
				stationaryIteration.Execute(rhsDense, solutionDense);
				++iter;

				convergenceMetric = convergenceCriterion.CalculateConvergenceMetric(solution, previousSolution);
				if (convergenceMetric < convergenceTolerance)
				{
					break;
				}
			}

			return new IterativeStatistics
			{
				AlgorithmName = stationaryIteration.Name,
				HasConverged = iter < maxIterations,
				NumIterationsRequired = iter,
				ConvergenceCriterion = (convergenceCriterion.DescribeConvergenceCriterion(convergenceTolerance), convergenceMetric),
			};
		}

		/// <summary>
		/// Constructs <see cref="StationaryIterativeMethod"/> instances, allows the user to specify some or all of the required 
		/// parameters and provides defaults for the rest.
		/// </summary>
		public class Factory
		{
			private readonly IStationaryIteration stationaryIteration;

			public Factory(IStationaryIteration stationaryIteration)
			{
				this.stationaryIteration = stationaryIteration;
			}

			public ISolutionConvergenceCriterion ConvergenceCriterion { get; set; } = new AbsoluteSolutionConvergenceCriterion();

			public double ConvergenceTolerance { get; set; } = 1E-10;


			/// <summary>
			/// Specifies how to calculate the maximum iterations that the GS algorithm will run for.
			/// </summary>
			public IMaxIterationsProvider MaxIterationsProvider { get; set; } = new PercentageMaxIterationsProvider(1.0);


			/// <summary>
			/// Creates a new instance of <see cref="StationaryIterativeMethod"/> with the specified parameters.
			/// </summary>
			public StationaryIterativeMethod Build()
			{
				return new StationaryIterativeMethod(stationaryIteration.CopyWithInitialSettings(), ConvergenceTolerance,
					ConvergenceCriterion.CopyWithInitialSettings(), MaxIterationsProvider.CopyWithInitialSettings());
			}
		}
	}
}
