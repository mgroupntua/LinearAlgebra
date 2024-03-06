using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;

using MGroup.LinearAlgebra.Commons;
using MGroup.LinearAlgebra.Distributed.Exceptions;
using MGroup.LinearAlgebra.Iterative;
using MGroup.LinearAlgebra.Iterative.ConjugateGradient;
using MGroup.LinearAlgebra.Iterative.Preconditioning;
using MGroup.LinearAlgebra.Iterative.Termination.Iterations;
using MGroup.LinearAlgebra.Matrices;
using MGroup.LinearAlgebra.Vectors;
using MGroup.MSolve.Solution.LinearSystem;

//TODO: Needs Builder pattern
//TODO: perhaps all quantities should be stored as mutable fields, exposed as readonly properties and the various strategies 
//      should read them from a reference of CG/PCG/PCPG, instead of having them injected.
//TODO: In regular CG, there is a check to perevent premature convergence, by correcting the residual. Can this be done for PCG 
//      as well? Would the preconditioned residual be updated as well?
namespace MGroup.LinearAlgebra.Distributed.IterativeMethods.PCG
{
	/// <summary>
	/// Implements the block Preconditioned Conjugate Gradient algorithm for solving linear systems with symmetric 
	/// positive definite matrices. This implementation is based on the algorithm presented in the s-step paper of Chronopoulos 
	/// </summary>
	public class BlockPcgAlgorithm : PcgAlgorithmBase
	{
		private const string name = "Block Preconditioned Conjugate Gradient";
		private readonly IPcgBetaParameterCalculation betaCalculation;
		private readonly IBlockPcgResidualUpdater blockPcgResidualUpdater;
		private BlockVectorOperator residualOperator;
		private int blockSize;

		private IGlobalVector[] residualKernels;
		private IGlobalVector[] directionKernels;
		private double[] residualSandwiches;	// 2n+1 sandwich products (r_i * M * r_j) of n-vector Krylov subspace (A*M, r)
		private double[] directionSandwiches;	// 2n+3 sandwich products (p_i * M * p_j) of n-vector Krylov subspace (A*M, p)
		private double[] residualDirectionSandwiches; // 2n+2 sandwich products (r_i * M * p_j) of vector Krylov subspace (A*M, r) with Krylov subspace (A*M, p)

		/// <summary>
		/// Block vector operator for calculating the residual of the current iteration
		/// </summary>
		internal BlockVectorOperator ResidualOperator { get => residualOperator; }
  
		/// <summary>
		/// Krylov subspace of (A * M)^n * r
		/// A is the matrix, M is inverse preconditioner matrix and r is the CG residual vector
		/// </summary>
		internal IGlobalVector[] ResidualKernels { get => residualKernels; }
  
		/// <summary>
		/// Krylov subspace of (A * M)^n * p
		/// A is the matrix, M is inverse preconditioner matrix and p is the CG conjugate direction vector
		/// </summary>
		internal IGlobalVector[] DirectionKernels { get => directionKernels; }

		private BlockPcgAlgorithm(int blockSize, double residualTolerance, IMaxIterationsProvider maxIterationsProvider,
			IPcgResidualConvergence pcgConvergence, IBlockPcgResidualUpdater residualUpdater, bool throwIfNotConvergence,
			IPcgBetaParameterCalculation betaCalculation) : 
			base(residualTolerance, maxIterationsProvider, pcgConvergence, residualUpdater, throwIfNotConvergence)
		{
			this.blockSize = blockSize;
			this.betaCalculation = betaCalculation;
			this.residualKernels = new IGlobalVector[blockSize];
			this.directionKernels = new IGlobalVector[blockSize + 1];
			this.residualSandwiches = new double[2 * blockSize - 1];
			this.directionSandwiches = new double[2 * blockSize + 1];
			this.residualDirectionSandwiches = new double[2 * blockSize];
			this.blockPcgResidualUpdater = residualUpdater;
		}

		/// <summary>
		/// Creates Krylov subspace kernel using the linear system matrix and its preconditioner, and the vector provided. 
		/// </summary>
		/// <param name="vector">
		/// The vector to be used for calculating the Krylov subspace kernel.
		/// </param>
		/// <param name="kernel">
		/// The array of vectors to which each kernel calculated will be stored at.
		/// </param>
		/// <remarks>
		/// Normally this function must be run in parallel in multiple kernels.
		/// </remarks>
		private void EvaluateKernel(IGlobalVector vector, IGlobalVector[] kernel)
		{
			kernel[0].CopyFrom(vector);
			for (int i = 1; i < kernel.Length; ++i)
			{
				var v1 = vector.CreateZero();
				Preconditioner.Apply(kernel[i - 1], v1);
				Matrix.MultiplyVector(v1, kernel[i]);
			}
		}

		/// <summary>
		/// Computes preconditioned dot products with the vectors or 2 Krylov subspaces.
		/// </summary>
		/// <param name="kernel1">The first kernel.</param>
		/// <param name="kernel2">The second kernel.</param>
		/// <param name="sandwich">The array to which the result is stored at.</param>
  		/// <remarks>
		/// Normally this function must be run in parallel in multiple kernels.
		/// If R(i) is the a vector of R Krylov subspace and P(i) a vector of P Krylov subspace,
		/// this function produces preconditioned dot products R(i) * M * P(i) where M is the inverse preconditioner matrix.
  		/// </remarks>
		private void EvaluateSandwich(IGlobalVector[] kernel1, IGlobalVector[] kernel2, double[] sandwich)
		{
			var v = kernel1[0].CreateZero();
			Preconditioner.Apply(kernel1[0], v);
			for (int i = 0; i < kernel2.Length; ++i)
				sandwich[i] = v.DotProduct(kernel2[i]);

			Preconditioner.Apply(kernel2[kernel2.Length - 1], v);
			for (int i = 1; i < kernel1.Length; ++i)
				sandwich[kernel2.Length + i - 1] = v.DotProduct(kernel1[i]);
		}

		/// <summary>
		/// Initializes information for the first block of algorithm.
		/// </summary>
  		/// <remarks>
		/// Normally this function must be run in parallel in multiple kernels.
  		/// </remarks>
		private void InitializeBlockInfo()
		{
			EvaluateKernel(residual, directionKernels);
			EvaluateSandwich(directionKernels, directionKernels, directionSandwiches);

			for (int i = 0; i < residualKernels.Length; ++i)    // Deep copy for vectors
				residualKernels[i].CopyFrom(directionKernels[i]);
			// shallow copy is enough for scalars
			Array.Copy(directionSandwiches, residualSandwiches, residualSandwiches.Length);
			Array.Copy(directionSandwiches, residualDirectionSandwiches, residualDirectionSandwiches.Length);
		}

		/// <summary>
		/// Updates information for the rest of the blocks of algorithm.
		/// </summary>
  		/// <remarks>
		/// Normally this function must be run in parallel in multiple kernels.
  		/// </remarks>
		private void UpdateBlockInfo()
		{
			EvaluateKernel(residual, residualKernels);
			EvaluateKernel(direction, directionKernels);
			EvaluateSandwich(residualKernels, residualKernels, residualSandwiches);
			EvaluateSandwich(directionKernels, directionKernels, directionSandwiches);
			EvaluateSandwich(residualKernels, directionKernels, residualDirectionSandwiches);
		}

		/// <summary>
		/// Evaluates the solution vector using the supplied block vector linear combination coefficients.
		/// </summary>
		/// <param name="solutionCoefficients">The block vector linear combination coefficients to be used for the calculation of the solution vector.</param>
  		private IGlobalVector EvaluateSolutionVector(BlockVectorOperator solutionCoefficients)
		{
			var x = residualKernels[0].CreateZero();
			Preconditioner.Apply(solutionCoefficients.EvaluateVector(residualKernels, directionKernels), x);
			return x;
		}

		protected override IterativeStatistics SolveInternal(int maxIterations, Func<IGlobalVector> zeroVectorInitializer)
		{
			//CalculateAndPrintExactResidual();

			// In contrast to the source algorithm, we initialize s here. At each iteration it will be overwritten, 
			// thus avoiding allocating & deallocating a new vector.
			//precondResidual = zeroVectorInitializer();

			for (int i = 0; i < residualKernels.Length; i++)
				residualKernels[i] = zeroVectorInitializer();
			for (int i = 0; i < directionKernels.Length; i++)
				directionKernels[i] = zeroVectorInitializer();

			InitializeBlockInfo();
			resDotPrecondRes = residualSandwiches[0]; // rr = r * M * r
			double initialResDotPrecondRes = resDotPrecondRes;

			// The convergence and beta strategies must be initialized immediately after the first r and r*inv(M)*r are computed.
			convergence.Initialize(this);
			betaCalculation.Initialize(this);

			// This is also used as output
			double residualNormRatio = double.NaN;

			for (iteration = 0; iteration < maxIterations; ++iteration)
			{
				var solutionOperator = new BlockVectorOperator(blockSize);
				residualOperator = new BlockVectorOperator(blockSize); 
				var directionOperator = new BlockVectorOperator(blockSize);

				residualOperator.R.Add(1d);
				directionOperator.P.Add(1d);

				for (var i = 0; i < blockSize; ++i)
				{
					double pAp = directionOperator.Sandwich(residualSandwiches, directionSandwiches, residualDirectionSandwiches); // pAp = p * M * A * M * p
					if (pAp <= 0)
					{
						if (i != 0)
						{
							solution.AddIntoThis(EvaluateSolutionVector(solutionOperator));  // x = x0 + M * x. It multiplied with M, because it should be
						}

						// We reached the max iterations before BlockPCG converged
						if (throwIfNotConvergence)
						{
							throw new IterativeMethodDidNotConvergeException("Block PCG terminated after the max allowable number of iterations =" +
								$" {maxIterations}, without reaching the required residual norm ratio tolerance = {ResidualTolerance}." +
								$" In contrast the final residual norm ratio is {residualNormRatio}. To accept solutions without convergence" +
								$" set BlockPcgAlgorithm.Builder.ThrowExceptionIfNotConvergence = false");
						}
						else
						{
							return new IterativeStatistics
							{
								HasConverged = false,
								HasStagnated = true,
								AlgorithmName = name,
								NumIterationsRequired = iteration * blockSize + i,
								ResidualNormRatioEstimation = resDotPrecondRes / initialResDotPrecondRes
							};
						}
					}

					stepSize = resDotPrecondRes / pAp;
					solutionOperator.UpdateX(stepSize, directionOperator);	// x += a * p, but in should be x += a * M * p
					residualOperator.UpdateR(stepSize, directionOperator);	// r -= a * A * p, but in should be r -= a * A * M * p

					resDotPrecondResOld = resDotPrecondRes;
					resDotPrecondRes = residualOperator.Square(residualSandwiches, directionSandwiches, residualDirectionSandwiches);	// rr2 = r * M * r

					// The default Fletcher-Reeves formula is: β = δnew / δold = (sNew * rNew) / (sOld * rOld)
					// However we could use a different one, e.g. for variable preconditioning Polak-Ribiere is usually better.
					paramBeta = betaCalculation.CalculateBeta(this);
					directionOperator.UpdateP(residualOperator, paramBeta); // p = r + b * p, but in should be p = r + b * M * p

					residualNormRatio = convergence.EstimateResidualNormRatio(this);
					if (residualNormRatio <= ResidualTolerance)
					{
						if (i != 0)
						{
							solution.AddIntoThis(EvaluateSolutionVector(solutionOperator));  // x = x0 + M * x. It multiplied with M, because it should be
						}

						return new IterativeStatistics
						{
							HasConverged = true,
							HasStagnated = false,
							AlgorithmName = name,
							NumIterationsRequired = iteration * blockSize + i,
							ResidualNormRatioEstimation = resDotPrecondRes / initialResDotPrecondRes,
						};
					}
				}

				// Evaluates solution (x).
				// x = x0 + M * x. It multiplied with M, because it should be:
				// on all previous calculations M was ommited
				solution.AddIntoThis(EvaluateSolutionVector(solutionOperator));
				// Evaluates direction (p).
				// It didn't multiplied with M, but it should be:
				// We use direction (p) and residual (r) non-multiplied with M inside block
				// This is the reason, why solution (x) needs multiplication with M after the block ends
				direction = directionOperator.EvaluateVector(residualKernels, directionKernels);
				// Evaluates residual (r).
				// It didn't multiplied with M, but it should be:
				// We use direction (p) and residual (r) non-multiplied with M inside block
				// This is the reason, why solution (x) needs multiplication with M after the block ends
				blockPcgResidualUpdater.UpdateResidual(this, residual);
				UpdateBlockInfo();
			}

			// We reached the max iterations before BlockPCG converged
			if (throwIfNotConvergence)
			{
				throw new IterativeMethodDidNotConvergeException("Block PCG terminated after the max allowable number of iterations =" +
					$" {maxIterations}, without reaching the required residual norm ratio tolerance = {ResidualTolerance}." +
					$" In contrast the final residual norm ratio is {residualNormRatio}. To accept solutions without convergence" +
					$" set BlockPcgAlgorithm.Builder.ThrowExceptionIfNotConvergence = false");
			}
			else
			{
				return new IterativeStatistics
				{
					AlgorithmName = name,
					HasConverged = false,
					NumIterationsRequired = maxIterations,
					ResidualNormRatioEstimation = residualNormRatio
				};
			}
		}

		private void CalculateAndPrintExactResidual()
		{
			var res = Rhs.CreateZero();
			Matrix.MultiplyVector(solution, res);
			res.SubtractIntoThis(Rhs);
			double norm = res.Norm2();
			Debug.WriteLine($"Iteration {iteration}: norm(r) = {norm}");
		}

		/// <summary>
		/// Constructs <see cref="BlockPcgAlgorithm"/> instances, allows the user to specify some or all of the required parameters 
		/// and provides defaults for the rest.
		/// </summary>
		public class Builder : PcgBuilderBase
		{
			/// <summary>
			/// Specifies how to calculate the beta parameter of PCG, which is used to update the direction vector. 
			/// </summary>
			public int BlockSize { get; set; } = 6;

			/// <summary>
			/// Specifies how often the residual vector will be corrected by an exact (but costly) calculation.
			/// This parameter returns <see cref="Builder.BlockResidualUpdater"/>
			/// </summary>
			public override IPcgResidualUpdater ResidualUpdater
			{
				get => BlockResidualUpdater;
				set
				{
					if (value is IBlockPcgResidualUpdater)
					{
						BlockResidualUpdater = (IBlockPcgResidualUpdater)value;
					}
					else
					{
						throw new ArgumentException($"Block PCG builder needs a descendent of {typeof(IBlockPcgResidualUpdater)} while value is {value.GetType()}.");
					}
				}
			}

			/// <summary>
			/// Specifies how often the residual vector will be corrected by an exact (but costly) calculation.
			/// </summary>
			public IBlockPcgResidualUpdater BlockResidualUpdater { get; set; } = new RegularBlockPcgResidualUpdater();

			/// <summary>
			/// Specifies how to calculate the beta parameter of PCG, which is used to update the direction vector. 
			/// </summary>
			public IPcgBetaParameterCalculation BetaCalculation { get; set; } = new FletcherReevesBeta();

			/// <summary>
			/// Creates a new instance of <see cref="BlockPcgAlgorithm"/>.
			/// </summary>
			public BlockPcgAlgorithm Build()
			{
				return new BlockPcgAlgorithm(BlockSize, ResidualTolerance, MaxIterationsProvider, Convergence,
											BlockResidualUpdater, ThrowExceptionIfNotConvergence, BetaCalculation);
			}
		}
	}
}
