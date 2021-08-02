using System;
using System.Diagnostics;
using MGroup.LinearAlgebra.Commons;
using MGroup.LinearAlgebra.Iterative.Preconditioning;
using MGroup.LinearAlgebra.Iterative.Termination;
using MGroup.LinearAlgebra.Vectors;

//TODO: I would rather implement reorthogonalization as an alternative strategy, rather than a different class.
//TODO: needs builder
namespace MGroup.LinearAlgebra.Iterative.PreconditionedConjugateGradient
{
	/// <summary>
	/// Implements the untransformed Preconditioned Conjugate Gradient algorithm for solving linear systems with symmetric 
	/// positive definite matrices. The implementation is based on the algorithm presented in pages 51-54 of the PhD dissertation 
	/// "Seismic soil-structure interaction with finite elements and the method of substructures", George Stavroulakis, 2014
	/// Authors: Serafeim Bakalakos, George Stavroulakis 
	/// </summary>
	public class ReorthogonalizedPcg : PcgAlgorithmBase
	{
		private const string name = "Reorthogonalized PCG";


		private ReorthogonalizedPcg(double residualTolerance, IMaxIterationsProvider maxIterationsProvider,
			IPcgResidualConvergence residualConvergence, IPcgResidualUpdater residualCorrection) :
			base(residualTolerance, maxIterationsProvider, residualConvergence, residualCorrection)
		{
			Convergence = residualConvergence; //TODO: Now there are 2 convergence properties. One here and one in base class. Fix it.
		}

		public IPcgResidualConvergence Convergence { get; set; }

		/// <summary>
		/// The dot product d * (A*d), where d is the direction vector <see cref="PcgAlgorithmBase.Direction"/>.
		/// </summary>
		public double DirectionTimesMatrixTimesDirection { get; private set; }

		//TODO: this could be abstracted to use a cyclic cache.
		public PcgReorthogonalizationCache ReorthoCache { get; set; } = new PcgReorthogonalizationCache();

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
			for (int i = 0; i < ReorthoCache.Directions.Count; ++i)
			{
				// x_d[i] = (d_i * b) / (d_i * q_i) 
				double xd = ReorthoCache.Directions[i].DotProduct(rhsNew) / ReorthoCache.DirectionsTimesMatrixTimesDirections[i];

				Debug.Assert(!double.IsNaN(xd));
				Debug.Assert(!double.IsPositiveInfinity(xd));
				Debug.Assert(!double.IsNegativeInfinity(xd));

				// x0 += d_i * x_d[i]
				initialSolution.AxpyIntoThis(ReorthoCache.Directions[i], xd);
			}
		}

		/// <summary>
		/// See <see cref="PcgAlgorithmBase.Clear"/>
		/// </summary>
		public override void Clear()
		{
			base.Clear();
			DirectionTimesMatrixTimesDirection = 0.0;
		}

		public override IterativeStatistics Solve(ILinearTransformation matrix, IPreconditioner preconditioner, IVectorView rhs,
			IVector solution, bool initialGuessIsZero, Func<IVector> zeroVectorInitializer)
		{
			//TODO: find a better way to handle optimizations for the case x0=0, than using an initialGuessIsZero flag
			Preconditions.CheckMultiplicationDimensions(matrix.NumColumns, solution.Length);
			Preconditions.CheckSystemSolutionDimensions(matrix.NumRows, rhs.Length);

			this.Matrix = matrix;
			this.Preconditioner = preconditioner;
			this.Rhs = rhs;

			// Initial solution and rhs (r = b - A * x)
			this.solution = solution;
			if (ReorthoCache.Directions.Count > 0)
			{
				if (!initialGuessIsZero) solution.Clear();
				CalculateInitialSolutionFromStoredDirections(rhs, solution);
				residual = ExactResidual.Calculate(matrix, rhs, solution);
			}
			else // preferably call base method
			{
				// r = b - A * x
				if (initialGuessIsZero) residual = rhs.Copy();
				else residual = ExactResidual.Calculate(matrix, rhs, solution);
			}

			// Initialize vectors 
			//TODO: Pehaps I can just clear them from previous iterations 
			precondResidual = zeroVectorInitializer();
			direction = zeroVectorInitializer();
			matrixTimesDirection = zeroVectorInitializer();

			int maxIterations = MaxIterationsProvider.GetMaxIterations(matrix.NumColumns);
			return SolveInternal(maxIterations, zeroVectorInitializer);
		}

		protected override IterativeStatistics SolveInternal(int maxIterations, Func<IVector> zeroVectorInitializer)
		{
			iteration = 0;
			Preconditioner.SolveLinearSystem(residual, precondResidual);

			// d0 = s0 = inv(M) * r0
			//direction.CopyFrom(precondResidual);
			//Preconditioner.SolveLinearSystem(residual, direction);
			UpdateDirectionVector(precondResidual, direction);

			// q0 = A * d0
			Matrix.Multiply(direction, matrixTimesDirection);
			DirectionTimesMatrixTimesDirection = direction.DotProduct(matrixTimesDirection);

			// Update the direction vectors cache
			ReorthoCache.StoreDirectionData(this);

			// δnew = δ0 = r0 * s0 = r0 * d0
			resDotPrecondRes = residual.DotProduct(direction);

			// The convergence strategy must be initialized immediately after the first r and r*inv(M)*r are computed.
			Convergence.Initialize(this);
			Stagnation.StoreInitialError(Convergence.EstimateResidualNormRatio(this));

			// This is also used as output
			double residualNormRatio = double.NaN;

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

				/// At this point we can check if CG has converged and exit, thus avoiding the uneccesary operations that follow.
				residualNormRatio = Convergence.EstimateResidualNormRatio(this);
				//Debug.WriteLine($"Reorthogonalized PCG iteration = {iteration}: residual norm ratio = {residualNormRatio}");
				Stagnation.StoreNewError(residualNormRatio);
				bool hasStagnated = Stagnation.HasStagnated();
				if (residualNormRatio <= ResidualTolerance)
				{
					return new IterativeStatistics
					{
						AlgorithmName = name,
						HasConverged = true,
						HasStagnated = false,
						NumIterationsRequired = iteration + 1,
						ResidualNormRatioEstimation = residualNormRatio
					};
				}
				if (hasStagnated)
				{
					return new IterativeStatistics
					{
						AlgorithmName = name,
						HasConverged = false,
						HasStagnated = true,
						NumIterationsRequired = iteration + 1,
						ResidualNormRatioEstimation = residualNormRatio
					};
				}

				// Update the direction vector using previous cached direction vectors.
				UpdateDirectionVector(precondResidual, direction);

				// q = A * d
				Matrix.Multiply(direction, matrixTimesDirection);
				DirectionTimesMatrixTimesDirection = direction.DotProduct(matrixTimesDirection);

				// Update the direction vectors cache
				ReorthoCache.StoreDirectionData(this);

				// α = (d * r) / (d * q) = (d * r) / (d * (A * d)) 
				stepSize = direction.DotProduct(residual) / DirectionTimesMatrixTimesDirection;
			}

			// We reached the max iterations before PCG converged
			return new IterativeStatistics
			{
				AlgorithmName = name,
				HasConverged = false,
				HasStagnated = false,
				NumIterationsRequired = maxIterations,
				ResidualNormRatioEstimation = residualNormRatio
			};
		}

		private void UpdateDirectionVector(IVectorView preconditionedResidual, IVector direction)
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
		}

		/// <summary>
		/// Constructs <see cref="ReorthogonalizedPcg"/> instances, allows the user to specify some or all of the 
		/// required parameters and provides defaults for the rest.
		/// Author: Serafeim Bakalakos
		/// </summary>
		public class Builder : PcgBuilderBase
		{
			public Builder()
			{
				Convergence = new RhsNormalizedConvergence();
			}

			/// <summary>
			/// Creates a new instance of <see cref="ReorthogonalizedPcg"/>.
			/// </summary>
			public ReorthogonalizedPcg Build()
			{
				return new ReorthogonalizedPcg(ResidualTolerance, MaxIterationsProvider, Convergence, ResidualUpdater);
			}
		}
	}
}
