//TODO: I would rather implement reorthogonalization as an alternative strategy, rather than a different class.
namespace MGroup.LinearAlgebra.Iterative.PreconditionedConjugateGradient.Reorthogonalization
{
	using System;
	using System.Diagnostics;

	using MGroup.LinearAlgebra.Commons;
	using MGroup.LinearAlgebra.Exceptions;
	using MGroup.LinearAlgebra.Iterative.Preconditioning;
	using MGroup.LinearAlgebra.Iterative.Termination.Iterations;
	using MGroup.LinearAlgebra.Iterative.Termination.Stegnation;
	using MGroup.LinearAlgebra.Vectors;

	/// <summary>
	/// Implements the untransformed Preconditioned Conjugate Gradient algorithm for solving linear systems with symmetric 
	/// positive definite matrices. The implementation is based on the algorithm presented in pages 51-54 of the PhD dissertation 
	/// "Seismic soil-structure interaction with finite elements and the method of substructures", George Stavroulakis, 2014
	/// </summary>
	public class ReorthogonalizedPcg : PcgAlgorithmBase
	{
		private const string name = "Reorthogonalized PCG";
		private readonly IPcgResidualUpdater residualUpdater;
		private readonly bool useDirectionVectorsOnlyForInitialSolution;

		private ReorthogonalizedPcg(double residualTolerance, IMaxIterationsProvider maxIterationsProvider,
			IPcgResidualConvergence residualConvergence, IPcgResidualUpdater residualUpdater, bool throwIfNotConvergence,
			IDirectionVectorsRetention directionVectorsRetention, bool useDirectionVectorsOnlyForInitialSolution)
			: base(residualTolerance, maxIterationsProvider, residualConvergence, throwIfNotConvergence)
		{
			this.residualUpdater = residualUpdater;
			this.useDirectionVectorsOnlyForInitialSolution = useDirectionVectorsOnlyForInitialSolution;
			DirectionVectorsRetention = directionVectorsRetention;
		}

		/// <summary>
		/// The dot product d * (A*d), where d is the direction vector <see cref="PcgAlgorithmBase.Direction"/>.
		/// </summary>
		public double DirectionTimesMatrixTimesDirection { get; private set; }

		public IDirectionVectorsRetention DirectionVectorsRetention { get; }

		//TODO: this could be abstracted to use a cyclic cache.
		public PcgReorthogonalizationCache ReorthoCache { get; set; } = new PcgReorthogonalizationCache();

		public double ResidualNormRatio { get; set; }

		public IStagnationCriterion Stagnation { get; set; } = new NullStagnationCriterion();

		/// <summary>
		/// Calculates the initial approximation to the linear system's solution vector, by using a series of conjugate direction 
		/// vectors that have been stored previously by PCG during the solution of other linear systems with the same matrix. 
		/// This method should be used to solve linear systems with different right hand side vectors and the same matrix.
		/// </summary>
		/// <param name="rhsNew">The right hand side vector of the new linear system.</param>
		/// <param name="initialSolution">
		/// The initial approximation to the solution vector, which PCG will improve. It will be overwritten by this method.
		/// </param>
		/// <exception cref="InvalidOperationException">Thrown if there are no direction vectors stored yet.</exception>
		public void CalculateInitialSolutionFromStoredDirections(IVectorView rhsNew, IVector initialSolution)
		{
			//TODO: An implementation by G. Stavroulakis discarded the last stored direction vector at this point. Why?
			//reorthoCache.RemoveNewDirectionVectorData(1);

			// x0 = D_nd * x_d, x_d = inv(Q_nd * D_nd) * D_nd^T * b
			// D_nd = [d_1 ... d_nd], Q_nd = A * D_nd = [q_1 ... q_nd], Q_nd * D_nd = diag([d1*A*d1 ... d_nd*A*d_nd])
			for (var i = 0; i < ReorthoCache.Directions.Count; ++i)
			{
				// x_d[i] = (d_i * b) / (d_i * q_i) 
				var xd = ReorthoCache.Directions[i].DotProduct(rhsNew) / ReorthoCache.DirectionsTimesMatrixTimesDirections[i];

				Debug.Assert(!double.IsNaN(xd));
				Debug.Assert(!double.IsPositiveInfinity(xd));
				Debug.Assert(!double.IsNegativeInfinity(xd));

				// x0 += d_i * x_d[i]
				initialSolution.AxpyIntoThis(ReorthoCache.Directions[i], xd);
			}
		}

		/// <summary>
		/// <inheritdoc/>
		/// </summary>
		public override void Clear()
		{
			base.Clear();
			DirectionTimesMatrixTimesDirection = 0.0;
			ReorthoCache.Clear();
		}

		public override IterativeStatistics Solve(ILinearTransformation matrix, IPreconditioner preconditioner, IVectorView rhs,
			IVector solution, bool initialGuessIsZero)
		{
			//TODO: find a better way to handle optimizations for the case x0=0, than using an initialGuessIsZero flag
			Preconditions.CheckMultiplicationDimensions(matrix.NumColumns, solution.Length);
			Preconditions.CheckSystemSolutionDimensions(matrix.NumRows, rhs.Length);

			Matrix = matrix;
			Preconditioner = preconditioner;
			Rhs = rhs;

			// Initial solution and rhs (r = b - A * x)
			this.solution = solution;
			if (ReorthoCache.Directions.Count > 0)
			{
				if (!initialGuessIsZero)
				{
					solution.Clear();
				}

				CalculateInitialSolutionFromStoredDirections(rhs, solution);
				if (useDirectionVectorsOnlyForInitialSolution)
				{
					ReorthoCache.Clear();
				}

				residual = ExactResidual.Calculate(matrix, rhs, solution);
			}
			else // preferably call base method
			{
				// r = b - A * x
				if (initialGuessIsZero)
				{
					residual = rhs.Copy();
				}
				else
				{
					residual = ExactResidual.Calculate(matrix, rhs, solution);
				}
			}

			// Initialize vectors
			//TODO: Pehaps I can just clear them from previous iterations 
			precondResidual = Rhs.CreateZeroVectorWithSameFormat();
			direction = Solution.CreateZeroVectorWithSameFormat();
			matrixTimesDirection = Rhs.CreateZeroVectorWithSameFormat();

			var maxIterations = MaxIterationsProvider.GetMaxIterations(matrix.NumColumns);
			ReorthoCache.StartGeneration();
			DirectionVectorsRetention.Initialize(this);

			IterativeStatistics stats = SolveInternal(maxIterations);

			DirectionVectorsRetention.DiscardDirectionVectors();
			return stats;
		}

		protected override IterativeStatistics SolveInternal(int maxIterations)
		{
			iteration = 0;

			// This is also used as output
			ResidualNormRatio = double.NaN;

			Preconditioner.SolveLinearSystem(residual, precondResidual);

			// Update the direction vector (d) and matrix-times-direction vector (q) and the stored vectors cache
			UpdateDirectionVector(precondResidual, direction);

			// δnew = δ0 = r0 * s0 = r0 * d0
			resDotPrecondRes = residual.DotProduct(direction);

			// The convergence strategy must be initialized immediately after the first r and r*inv(M)*r are computed.
			convergence.Initialize(this);
			Stagnation.StoreInitialError(convergence.EstimateResidualNormRatio(this));

			// α0 = (d0 * r0) / (d0 * q0) = (s0 * r0) / (d0 * (A * d0))
			stepSize = resDotPrecondRes / DirectionTimesMatrixTimesDirection;

			for (iteration = 1; iteration < maxIterations; ++iteration)
			{
				// x = x + α * d
				solution.AxpyIntoThis(direction, stepSize);

				// Normally the residual vector is updated as: r = r - α * q. However corrections might need to be applied.
				residualUpdater.UpdateResidual(this, residual);

				// s = inv(M) * r
				Preconditioner.SolveLinearSystem(residual, precondResidual);

				// δold = δnew
				resDotPrecondResOld = resDotPrecondRes;

				// δnew = r * s
				resDotPrecondRes = residual.DotProduct(precondResidual);

				// At this point we can check if CG has converged and exit, thus avoiding the uneccesary operations that follow.
				ResidualNormRatio = convergence.EstimateResidualNormRatio(this);
				//Debug.WriteLine($"Reorthogonalized PCG iteration = {iteration}: residual norm ratio = {residualNormRatio}");
				Stagnation.StoreNewError(ResidualNormRatio);
				bool hasStagnated = Stagnation.HasStagnated();
				if (hasStagnated || (ResidualNormRatio <= ResidualTolerance))
				{
					return new IterativeStatistics
					{
						AlgorithmName = name,
						HasConverged = true,
						HasStagnated = hasStagnated,
						NumIterationsRequired = iteration + 1,
						ResidualNormRatioEstimation = ResidualNormRatio,
					};
				}

				// Update the direction vector (d) and matrix-times-direction vector (q) and the stored vectors cache
				UpdateDirectionVector(precondResidual, direction);

				// α = (d * r) / (d * q) = (d * r) / (d * (A * d)) 
				stepSize = direction.DotProduct(residual) / DirectionTimesMatrixTimesDirection;
			}

			// We reached the max iterations before PCG converged
			if (throwIfNotConvergence)
			{
				throw new IterativeMethodDidNotConvergeException(
					"Reorthogonalized-PCG terminated after the max allowable number of iterations =" +
					$" {maxIterations}, without reaching the required residual norm ratio tolerance = {ResidualTolerance}." +
					$" In contrast the final residual norm ratio is {ResidualNormRatio}. To accept solutions without convergence" +
					$" set ReorthogonalizedPcg.Factory.ThrowExceptionIfNotConvergence = false");
			}
			else
			{
				return new IterativeStatistics
				{
					AlgorithmName = name,
					HasConverged = false,
					HasStagnated = false,
					NumIterationsRequired = maxIterations,
					ResidualNormRatioEstimation = ResidualNormRatio,
				};
			}
		}

		private void UpdateDirectionVector(IVectorView preconditionedResidual, IVector direction)
		{
			bool useReortho = DirectionVectorsRetention.KeepUsingReorthogonalization();

			if (useReortho)
			{
				// d = s - sum(β_i * d_i), 0 <= i < numStoredDirections
				// β_i = (s * q_i) / (d_i * q_i)
				direction.CopyFrom(preconditionedResidual);
				for (int i = 0; i < ReorthoCache.Directions.Count; ++i)
				{
					double beta = preconditionedResidual.DotProduct(ReorthoCache.MatrixTimesDirections[i])
						/ ReorthoCache.DirectionsTimesMatrixTimesDirections[i];
					direction.AxpyIntoThis(ReorthoCache.Directions[i], -beta);
				}

				// q = A * d
				Matrix.Multiply(direction, matrixTimesDirection);
				DirectionTimesMatrixTimesDirection = direction.DotProduct(matrixTimesDirection);

				// Update the direction vectors cache
				ReorthoCache.StoreDirectionData(this);
			}
			else
			{
				// Direction vector update based on Fletcher-Reeves formula, like in regular PCG
				// β = δnew / δold = (sNew * rNew) / (sOld * rOld)
				paramBeta = resDotPrecondRes / resDotPrecondResOld;

				// d = s + β * d
				direction.LinearCombinationIntoThis(paramBeta, precondResidual, 1.0);

				// q = A * d
				Matrix.Multiply(direction, matrixTimesDirection);
				DirectionTimesMatrixTimesDirection = direction.DotProduct(matrixTimesDirection);
			}
		}

		/// <summary>
		/// Constructs <see cref="ReorthogonalizedPcg"/> instances, allows the user to specify some or all of the 
		/// required parameters and provides defaults for the rest.
		/// </summary>
		public class Factory : PcgFactoryBase
		{
			public Factory()
			{
				Convergence = new PureResidualConvergence();
				DirectionVectorsRetention = new PercentageDirectionVectorsRetention(1.1);
				UseDirectionVectorsOnlyForInitialSolution = false;
			}

			public IDirectionVectorsRetention DirectionVectorsRetention { get; set; }

			/// <summary>
			/// Specifies how often the residual vector will be corrected by an exact (but costly) calculation.
			/// </summary>
			public virtual IPcgResidualUpdater ResidualUpdater { get; set; } = new DefaultPcgResidualUpdater();

			public bool UseDirectionVectorsOnlyForInitialSolution { get; set; }

			/// <summary>
			/// Creates a new instance of <see cref="ReorthogonalizedPcg"/>.
			/// </summary>
			public ReorthogonalizedPcg Build()
			{
				return new ReorthogonalizedPcg(ResidualTolerance, MaxIterationsProvider.CopyWithInitialSettings(),
					Convergence.CopyWithInitialSettings(), ResidualUpdater.CopyWithInitialSettings(),
					ThrowExceptionIfNotConvergence, DirectionVectorsRetention.CopyWithInitialSettings(),
					UseDirectionVectorsOnlyForInitialSolution);
			}
		}
	}
}
