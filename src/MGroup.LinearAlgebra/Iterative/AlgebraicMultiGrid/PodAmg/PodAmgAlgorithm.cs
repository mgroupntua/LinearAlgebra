namespace MGroup.LinearAlgebra.Iterative.AlgebraicMultiGrid.PodAmg
{
	using System;
	using System.Collections.Generic;
	using System.Text;
	using System.Xml.Linq;

	using MGroup.LinearAlgebra.Commons;
	using MGroup.LinearAlgebra.Iterative.AlgebraicMultiGrid.Smoothing;
	using MGroup.LinearAlgebra.Iterative.GaussSeidel;
	using MGroup.LinearAlgebra.Iterative.Preconditioning;
	using MGroup.LinearAlgebra.Iterative.Termination.Convergence;
	using MGroup.LinearAlgebra.Iterative.Termination.Iterations;
	using MGroup.LinearAlgebra.Matrices;
	using MGroup.LinearAlgebra.Triangulation;
	using MGroup.LinearAlgebra.Vectors;

	public class PodAmgAlgorithm : IDisposable
	{
		private const string name = "POD-AMG";

		private readonly IMaxIterationsProvider _maxIterationsProvider;
		private readonly ISolutionConvergenceCriterion _convergenceCriterion;
		private readonly double _convergenceTolerance;

		private readonly IMultigridSmoother _smoother;

		//TODO: encapsulate these. They are POD inputs.
		private readonly int _numPrincipalComponents;
		private readonly Matrix _sampleVectors;

		private CsrMatrix _fineMatrix;
		private CholeskyFull _coarseMatrixFactorized;

		// Restriction is the transpose of this
		private Matrix _prolongation; 

		private PodAmgAlgorithm(Matrix sampleVectors, int numPrincipalComponents, IMultigridSmoother smoother,
			double convergenceTolerance, ISolutionConvergenceCriterion convergenceCriterion, 
			IMaxIterationsProvider maxIterationsProvider) 
		{
			_sampleVectors = sampleVectors;
			_numPrincipalComponents = numPrincipalComponents;
			_smoother = smoother;
			this._convergenceTolerance = convergenceTolerance;
			_convergenceCriterion = convergenceCriterion;
			_maxIterationsProvider = maxIterationsProvider;
		}

		public void Dispose() 
		{ 
			if (_smoother != null) _smoother.Dispose();
		}

		public void Initialize(CsrMatrix systemMatrix)
		{
			if (_prolongation == null) // It may have been created in previous system solutions
			{
				var pod = new ProperOrthogonalDecomposition();
				_prolongation = pod.CalculatePrincipalComponents(_sampleVectors.NumColumns, _sampleVectors, _numPrincipalComponents);
			}

			_fineMatrix = systemMatrix;
			_smoother.Initialize(_fineMatrix);

			Matrix temp = _fineMatrix.MultiplyRight(_prolongation);
			Matrix coarseMatrix = _prolongation.MultiplyRight(temp, transposeThis:true, transposeOther:false);

			_coarseMatrixFactorized = CholeskyFull.Factorize(coarseMatrix.NumRows, coarseMatrix.RawData);
		}

		/// <summary>
		/// </summary>
		/// <param name="rhs"></param>
		/// <param name="solution">An initial guess or a zero vector.</param>
		/// <returns></returns>
		public IterativeStatistics Solve(Vector rhs, Vector solution)
		{
			Preconditions.CheckSquareLinearSystemDimensions(_fineMatrix, rhs, solution);
			int n0 = _fineMatrix.NumRows;
			var r0 = Vector.CreateZero(n0);
			var e0 = Vector.CreateZero(n0);
			int n1 = _coarseMatrixFactorized.Order;
			var r1 = Vector.CreateZero(n1);
			var e1 = Vector.CreateZero(n1);
			Vector previousX0 = solution.Copy();

			// Determine termination criteria
			double relativeTolerance = _convergenceTolerance * rhs.Norm2();
			relativeTolerance = (relativeTolerance == 0.0) ? _convergenceTolerance : relativeTolerance;
			double convergenceMetric = double.NaN;
			int maxCycles = _maxIterationsProvider.GetMaxIterations(n0);
			int iter = 0;
			
			while (iter < maxCycles)
			{
				previousX0.CopyFrom(solution);

				// Pre-smoothing on lvl 0 to get an estimate of the solution x0. Use the x0 from previous cycles as initial guess.
				_smoother.Smooth(rhs, solution);

				// Find the residual on lvl 0: r0=b-A0*x0
				//TODO: Use ExactResidual class for this
				_fineMatrix.MultiplyIntoResult(solution, r0);
				r0.LinearCombinationIntoThis(-1.0, rhs, 1.0);

				// Restrict lvl 0 residual to lvl 1: r1 = P^T * r0
				_prolongation.MultiplyIntoResult(r0, r1, transposeThis: true);

				// Find an estimate of the error on lvl 1 by solving exactly the system: A1*e1=r1.
				_coarseMatrixFactorized.SolveLinearSystem(r1, e1);

				// Interpolate the lvl 1 error estimate to lvl 0: e0 = P * e1
				_prolongation.MultiplyIntoResult(e1, e0, transposeThis: false);

				// Correct the solution estimate on lvl 0 using the interpolated error: x0 = x0 + e0
				solution.AddIntoThis(e0);

				// Post-smoothing on lvl 0 to further improve the solution estimate x0. Use the corrected x0 as initial guess.
				_smoother.Smooth(rhs, solution);

				// Check termination
				++iter;
				_fineMatrix.MultiplyIntoResult(solution, r0);
				r0.LinearCombinationIntoThis(-1.0, rhs, 1.0);
				convergenceMetric = r0.Norm2();
                System.Diagnostics.Debug.WriteLine($"Iter {iter-1}: norm2(r) = {convergenceMetric}");

				//convergenceMetric = _convergenceCriterion.CalculateConvergenceMetric(solution, previousX0);
				//if (convergenceMetric < _convergenceTolerance)
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
				ConvergenceCriterion = ("norm2(b - A * x) / norm2(b) < " + _convergenceTolerance, convergenceMetric)
				//ConvergenceCriterion =
				//	(_convergenceCriterion.DescribeConvergenceCriterion(_convergenceTolerance), convergenceMetric),
			};
		}

		public class Builder
		{
			public Builder()
			{
				// Defaults are taken from PyAMG
				ConvergenceCriterion = new AbsoluteSolutionConvergenceCriterion();
				ConvergenceTolerance = 1E-5;
				MaxIterationsProvider = new FixedMaxIterationsProvider(100);
				SmootherBuilder = new GaussSeidelSmoother.Builder(
					new GaussSeidelIterationCsrSerial.Builder(), 
					GaussSeidelSweepDirection.Forward, 
					numIterations: 1);
			}

			public ISolutionConvergenceCriterion ConvergenceCriterion { get; set; } 

			public double ConvergenceTolerance { get; set; }

			/// <summary>
			/// Specifies how to calculate the maximum iterations that the algorithm will run for.
			/// </summary>
			public IMaxIterationsProvider MaxIterationsProvider { get; set; }

			public IMultigridSmootherBuilder SmootherBuilder { get; set; }

			public PodAmgAlgorithm Create(Matrix sampleVectors, int numPrincipalComponents)
			{
				return new PodAmgAlgorithm(sampleVectors, numPrincipalComponents, SmootherBuilder.Create(), 
					ConvergenceTolerance, ConvergenceCriterion, MaxIterationsProvider);
			}
		}

	}
}
