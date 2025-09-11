namespace MGroup.LinearAlgebra.AlgebraicMultiGrid.PodAmg
{
	using MGroup.LinearAlgebra.AlgebraicMultiGrid;
	using MGroup.LinearAlgebra.Commons;
	using MGroup.LinearAlgebra.Iterative;
	using MGroup.LinearAlgebra.Iterative.Stationary.CSR;
	using MGroup.LinearAlgebra.Iterative.Termination.Convergence;
	using MGroup.LinearAlgebra.Iterative.Termination.Iterations;
	using MGroup.LinearAlgebra.Matrices;
	using MGroup.LinearAlgebra.Triangulation;
	using MGroup.LinearAlgebra.Vectors;

	/// <summary>
	/// Solves linear systems using the 2-level POD-AMG (Proper Orthogonal Decomposition - Algebraic Multigrid) method.
	/// </summary>
	public class PodAmgAlgorithm
	{
		private const string name = "POD-AMG";

		private readonly IMaxIterationsProvider maxIterationsProvider;
		private readonly ISolutionConvergenceCriterion convergenceCriterion;
		private readonly double convergenceTolerance;

		private readonly MultigridLevelSmoothing smoothing;

		//TODO: encapsulate these. They are POD inputs.
		private readonly bool keepOnlyNonZeroPrincipalComponents;
		private readonly int numPrincipalComponents;
		private readonly Matrix sampleVectors;

		private CsrMatrix fineMatrix;
		private CholeskyFull coarseMatrixFactorized;

		// Restriction is the transpose of this
		private Matrix _prolongation;

		private PodAmgAlgorithm(Matrix sampleVectors, bool keepOnlyNonZeroPrincipalComponents, int numPrincipalComponents,
			MultigridLevelSmoothing smoothing, double convergenceTolerance, ISolutionConvergenceCriterion convergenceCriterion,
			IMaxIterationsProvider maxIterationsProvider)
		{
			this.keepOnlyNonZeroPrincipalComponents = keepOnlyNonZeroPrincipalComponents;
			this.sampleVectors = sampleVectors;
			this.numPrincipalComponents = numPrincipalComponents;
			this.smoothing = smoothing;
			this.convergenceTolerance = convergenceTolerance;
			this.convergenceCriterion = convergenceCriterion;
			this.maxIterationsProvider = maxIterationsProvider;
		}

		/// <summary>
		/// Performs POD and prepares the multigrid operators for the provided matrix and training data.
		/// </summary>
		/// <param name="systemMatrix">The matrix of the linear system in CSR format.</param>
		public void Initialize(CsrMatrix systemMatrix)
		{
			if (_prolongation == null) // It may have been created in previous system solutions
			{
				var pod = new ProperOrthogonalDecomposition(keepOnlyNonZeroPrincipalComponents);
				_prolongation = pod.CalculatePrincipalComponents(sampleVectors.NumColumns, sampleVectors, numPrincipalComponents);
			}

			fineMatrix = systemMatrix;
			smoothing.UpdateMatrix(fineMatrix, true);

			var temp = fineMatrix.MultiplyRight(_prolongation);
			var coarseMatrix = _prolongation.MultiplyRight(temp, transposeThis: true, transposeOther: false);

			coarseMatrixFactorized = CholeskyFull.Factorize(coarseMatrix.NumRows, coarseMatrix.RawData);
		}

		/// <summary>
		/// Solves the linear system using the 2-level POD-AMG method. <see cref="Initialize(CsrMatrix)"/> must be called first.
		/// </summary>
		/// <param name="rhs">The right-hand-side vector of the linear system</param>
		/// <param name="solution">Initial guess for the solution vector. Can be zero.</param>
		/// <returns>Info about the performance of the algorithm.</returns>
		public IterativeStatistics Solve(Vector rhs, Vector solution)
		{
			Preconditions.CheckSquareLinearSystemDimensions(fineMatrix, rhs, solution);
			var n0 = fineMatrix.NumRows;
			var r0 = Vector.CreateZero(n0);
			var e0 = Vector.CreateZero(n0);
			var n1 = coarseMatrixFactorized.Order;
			var r1 = Vector.CreateZero(n1);
			var e1 = Vector.CreateZero(n1);
			var previousX0 = solution.Copy();

			// Determine termination criteria
			var relativeTolerance = convergenceTolerance * rhs.Norm2();
			relativeTolerance = relativeTolerance == 0.0 ? convergenceTolerance : relativeTolerance;
			var convergenceMetric = double.NaN;
			var maxCycles = maxIterationsProvider.GetMaxIterations(n0);
			var iter = 0;

			while (iter < maxCycles)
			{
				previousX0.CopyFrom(solution);

				// Pre-smoothing on lvl 0 to get an estimate of the solution x0. Use the x0 from previous cycles as initial guess.
				smoothing.ApplyPreSmoothers(rhs, solution);

				// Find the residual on lvl 0: r0=b-A0*x0
				//TODO: Use ExactResidual class for this
				fineMatrix.MultiplyIntoResult(solution, r0);
				r0.LinearCombinationIntoThis(-1.0, rhs, 1.0);

				// Restrict lvl 0 residual to lvl 1: r1 = P^T * r0
				_prolongation.MultiplyIntoResult(r0, r1, transposeThis: true);

				// Find an estimate of the error on lvl 1 by solving exactly the system: A1*e1=r1.
				coarseMatrixFactorized.SolveLinearSystem(r1, e1);

				// Interpolate the lvl 1 error estimate to lvl 0: e0 = P * e1
				_prolongation.MultiplyIntoResult(e1, e0, transposeThis: false);

				// Correct the solution estimate on lvl 0 using the interpolated error: x0 = x0 + e0
				solution.AddIntoThis(e0);

				// Post-smoothing on lvl 0 to further improve the solution estimate x0. Use the corrected x0 as initial guess.
				smoothing.ApplyPostSmoothers(rhs, solution);

				// Check termination
				++iter;
				fineMatrix.MultiplyIntoResult(solution, r0);
				r0.LinearCombinationIntoThis(-1.0, rhs, 1.0);
				convergenceMetric = r0.Norm2();
				System.Diagnostics.Debug.WriteLine($"Iter {iter - 1}: norm2(r) = {convergenceMetric}");

				//convergenceMetric = convergenceCriterion.CalculateConvergenceMetric(solution, previousX0);
				//if (convergenceMetric < convergenceTolerance)
				if (convergenceMetric <= relativeTolerance)
				{
					break;
				}
			}

			return new IterativeStatistics()
			{
				AlgorithmName = name,
				HasConverged = iter < maxCycles,
				NumIterationsRequired = iter,
				ConvergenceCriterion = ("norm2(b - A * x) / norm2(b) < " + convergenceTolerance, convergenceMetric)
				//ConvergenceCriterion =
				//	(convergenceCriterion.DescribeConvergenceCriterion(convergenceTolerance), convergenceMetric),
			};
		}

		/// <summary>
		/// Helper class to create instances of <see cref="PodAmgAlgorithm"/>.
		/// </summary>
		public class Factory
		{
			/// <summary>
			/// Creates new instance of <see cref="Factory"/>.
			/// </summary>
			public Factory()
			{
				KeepOnlyNonZeroPrincipalComponents = true;

				// Defaults are taken from PyAMG
				ConvergenceCriterion = new AbsoluteSolutionConvergenceCriterion();
				ConvergenceTolerance = 1E-5;
				MaxIterationsProvider = new FixedMaxIterationsProvider(100);
				Smoothing = new MultigridLevelSmoothing()
					.AddPreSmoother(new GaussSeidelIterationCsr(), 1)
					.AddPostSmoother(new GaussSeidelIterationCsr(), 1);
			}

			/// <summary>
			/// Specifies how to determines when the iterative algorithm has converged to the desired solution.
			/// </summary>
			public ISolutionConvergenceCriterion ConvergenceCriterion { get; set; }

			/// <summary>
			/// Max absolute value for the norm-2 of the residual vector to determine that the algorithm has converged.
			/// </summary>
			public double ConvergenceTolerance { get; set; }

			/// <summary>
			/// If true, eigenvalues that are smaller than a specified tolerance will be ignored during POD,
			/// along with their corresponding eigenvectors.
			/// </summary>
			public bool KeepOnlyNonZeroPrincipalComponents { get; set; }

			/// <summary>
			/// Specifies how to calculate the maximum iterations that the algorithm will run for.
			/// </summary>
			public IMaxIterationsProvider MaxIterationsProvider { get; set; }

			/// <summary>
			/// Specifies the smoothing operator (e.g. Gauss-Seidel, Jacobi, SOR, ...) of the multigrid procedure.
			/// </summary>
			public MultigridLevelSmoothing Smoothing { get; set; }

			/// <summary>
			/// Creates a new instance of <see cref="PodAmgAlgorithm"/> with the specified settings.
			/// </summary>
			/// <param name="sampleVectors">Matrix whose columns are the vectors to be used in POD.</param>
			/// <param name="numPrincipalComponents">How many principal components to keep during POD.</param>
			/// <returns>A new instance of <see cref="PodAmgAlgorithm"/>.</returns>
			public PodAmgAlgorithm Create(Matrix sampleVectors, int numPrincipalComponents)
			{
				return new PodAmgAlgorithm(sampleVectors, KeepOnlyNonZeroPrincipalComponents, numPrincipalComponents,
					Smoothing.CopyWithInitialSettings(), ConvergenceTolerance, ConvergenceCriterion.CopyWithInitialSettings(),
					MaxIterationsProvider.CopyWithInitialSettings());
			}
		}
	}
}
