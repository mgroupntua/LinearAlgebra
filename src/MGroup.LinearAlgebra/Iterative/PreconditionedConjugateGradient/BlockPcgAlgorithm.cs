using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;

using DotNumerics.Optimization.TN;

using MGroup.LinearAlgebra.Commons;
using MGroup.LinearAlgebra.Iterative.ConjugateGradient;
using MGroup.LinearAlgebra.Iterative.Preconditioning;
using MGroup.LinearAlgebra.Iterative.Termination.Iterations;
using MGroup.LinearAlgebra.Matrices;
using MGroup.LinearAlgebra.Vectors;

//TODO: Needs Builder pattern
//TODO: perhaps all quantities should be stored as mutable fields, exposed as readonly properties and the various strategies 
//      should read them from a reference of CG/PCG/PCPG, instead of having them injected.
//TODO: In regular CG, there is a check to perevent premature convergence, by correcting the residual. Can this be done for PCG 
//      as well? Would the preconditioned residual be updated as well?
namespace MGroup.LinearAlgebra.Iterative.PreconditionedConjugateGradient
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

		private IVector[] residualKernels;
		private IVector[] directionKernels;
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
		internal IVector[] ResidualKernels { get => residualKernels; }
		/// <summary>
		/// Krylov subspace of (A * M)^n * p
		/// A is the matrix, M is inverse preconditioner matrix and p is the CG conjugate direction vector
		/// </summary>
		internal IVector[] DirectionKernels { get => directionKernels; }

		private BlockPcgAlgorithm(int blockSize, double residualTolerance, IMaxIterationsProvider maxIterationsProvider,
			IPcgResidualConvergence pcgConvergence, IBlockPcgResidualUpdater residualUpdater, 
			IPcgBetaParameterCalculation betaCalculation) : 
			base(residualTolerance, maxIterationsProvider, pcgConvergence, residualUpdater)
		{
			this.blockSize = blockSize;
			this.betaCalculation = betaCalculation;
			this.residualKernels = new IVector[blockSize];
			this.directionKernels = new IVector[blockSize + 1];
			this.residualSandwiches = new double[2 * blockSize - 1];
			this.directionSandwiches = new double[2 * blockSize + 1];
			this.residualDirectionSandwiches = new double[2 * blockSize];
			this.blockPcgResidualUpdater = residualUpdater;
		}

		/// <summary>
		/// Creates Krylov subspace kernel using the linear system matrix and its preconditioner, and the vector provided. 
		/// Normally this function must be run in parallel in multiple kernels.
		/// </summary>
		/// <param name="vector">
		/// The vector to be used for calculating the Krylov subspace kernel.
		/// </param>
		/// <param name="kernel">
		/// The array of vectors to which each kernel calculated will be stored at.
		/// </param>
		private void EvaluateKernel(IVector vector, IVector[] kernel)
		{
			kernel[0].CopyFrom(vector);
			for (int i = 1; i < kernel.Length; ++i)
			{
				var v1 = vector.CreateZeroVectorWithSameFormat();
				Preconditioner.SolveLinearSystem(kernel[i - 1], v1);
				Matrix.Multiply(v1, kernel[i]);
			}
		}

		/// <summary>
		/// Creates preconditioned dot products with the vectors or 2 Krylov subspaces. 
		/// Normally this function must be run in parallel in multiple kernels.
		/// </summary>
		/// <param name="vector">
		/// The vector to be used for calculating the Krylov subspace kernel.
		/// </param>
		/// <param name="kernel">
		/// The array of vectors to which each kernel calculated will be stored at.
		/// </param>
		/// <param name="size">
		/// The power to which the kernel will be calculated
		/// </param>

		/// Creates preconditioned dot products with the vectors or 2 Krylov subspaces.
		/// Normally this function must be run in parallel in multiple kernels.
		/// If R(i) is the a vector of R Krylov subspace and P(i) a vector of P Krylov subspace,
		/// this function produces preconditioned dot products R(i) * M * P(i) where M is the inverse preconditioner matrix.
		/// <param name="preconditioner">The preconditioner.</param>
		/// <param name="krylov1">The krylov subspace R. A list of n vectors. First vector is vector (A * M)^0 * r
		/// and last vector is (A * M)^(n-1) * r, where A is the matrix and M is the inverse preconditioner.</param>.
		/// <param name="krylov2">The krylov subspace P. A list of n vectors. First vector is vector (A * M)^0 * p
		/// and last vector is (A * M)^(n-1) * p, where A is the matrix and M is the inverse preconditioner.</param>.
		/// <return>The result vector with the preconditioned dot products. This is the output of this function.</return>
		private void EvaluateSandwich(IVector[] kernel1, IVector[] kernel2, double[] sandwich)
		{
			var v = kernel1[0].CreateZeroVectorWithSameFormat();
			Preconditioner.SolveLinearSystem(kernel1[0], v);
			for (int i = 0; i < kernel2.Length; ++i)
				sandwich[i] = v.DotProduct(kernel2[i]);

			Preconditioner.SolveLinearSystem(kernel2[kernel2.Length - 1], v);
			for (int i = 1; i < kernel1.Length; ++i)
				sandwich[kernel2.Length + i - 1] = v.DotProduct(kernel1[i]);
		}

		/** Create all the information for the first block of algorithm.
 		 * Normally this function must be run in parallel in multiple kernels.
 		 * <param name="linearTransformation">A linear transformation (A * M) * vector where A is the matrix and M is the inverse preconditioner matrix.</param>
 		 * <param name="preconditioner">The preconditioner.</param>
 		 * <param name="r">Vector r of CG.</param>
		 * <param name="n">Number of vectors in Krylov subspace.</param>
 		 * <return>The results of parallel processing.</return>
 		 */
		private void UpdateBlockInfoUsingResidualOnly()
		{
			EvaluateKernel(residual, directionKernels);
			EvaluateSandwich(directionKernels, directionKernels, directionSandwiches);

			for (int i = 0; i < residualKernels.Length; ++i)    // Deep copy for vectors
				residualKernels[i].CopyFrom(directionKernels[i]);
			// shallow copy is enough for scalars
			Array.Copy(directionSandwiches, residualSandwiches, residualSandwiches.Length);
			Array.Copy(directionSandwiches, residualDirectionSandwiches, residualDirectionSandwiches.Length);
		}

		/** Create all the information for all the blocks of algorithm except first.
 		 * Normally this function must be run in parallel in multiple kernels.
 		 * <param name="linearTransformation">A linear transformation (A * M) * vector where A is the matrix and M is the inverse preconditioner matrix.</param>
 		 * <param name="preconditioner">The preconditioner.</param>
 		 * <param name="r">Vector r of CG.</param>
 		 * <param name="p">Vector p of CG.</param>
 		 * <param name="n">Number of vectors in Krylov subspace.</param>
 		 * <return>The results of parallel processing.</return>
		 */
		private void UpdateBlockInfo()
		{
			EvaluateKernel(residual, residualKernels);
			EvaluateKernel(direction, directionKernels);
			EvaluateSandwich(residualKernels, residualKernels, residualSandwiches);
			EvaluateSandwich(directionKernels, directionKernels, directionSandwiches);
			EvaluateSandwich(residualKernels, directionKernels, residualDirectionSandwiches);
		}

		private IVector EvaluateSolutionVector(BlockVectorOperator solutionCoefficients)
		{
			var x = residualKernels[0].CreateZeroVectorWithSameFormat();
			Preconditioner.SolveLinearSystem(solutionCoefficients.EvaluateVector(residualKernels, directionKernels), x);
			return x;
		}

		protected override IterativeStatistics SolveInternal(int maxIterations, Func<IVector> zeroVectorInitializer)
		{
			//CalculateAndPrintExactResidual();

			// In contrast to the source algorithm, we initialize s here. At each iteration it will be overwritten, 
			// thus avoiding allocating & deallocating a new vector.
			//precondResidual = zeroVectorInitializer();

			for (int i = 0; i < residualKernels.Length; i++)
				residualKernels[i] = zeroVectorInitializer();
			for (int i = 0; i < directionKernels.Length; i++)
				directionKernels[i] = zeroVectorInitializer();

			UpdateBlockInfoUsingResidualOnly();
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

						return new BlockIterativeStatistics
						{
							HasConverged = false,
							HasStagnated = true,
							AlgorithmName = name,
							NumIterationsRequired = iteration * blockSize + i,
							ResidualNormRatioEstimation = resDotPrecondRes / initialResDotPrecondRes,
							ResDotRes = resDotPrecondRes,
							Solution = solution,
						};
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

						return new BlockIterativeStatistics
						{
							HasConverged = true,
							HasStagnated = false,
							AlgorithmName = name,
							NumIterationsRequired = iteration * blockSize + i,
							ResidualNormRatioEstimation = resDotPrecondRes / initialResDotPrecondRes,
							ResDotRes = resDotPrecondRes,
							Solution = solution,
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
				// Because it updates residual which is by reference the residualKernel[0],
				// we run this after evaluation of solution (x) and direction (p),
				// otherwise evaluation will use destroyed residualKernel.
				// It didn't multiplied with M, but it should be:
				// We use direction (p) and residual (r) non-multiplied with M inside block
				// This is the reason, why solution (x) needs multiplication with M after the block ends
				blockPcgResidualUpdater.UpdateResidual(this, residual);
				UpdateBlockInfo();
			}

			// We reached the max iterations before converging
			return new IterativeStatistics
			{
				AlgorithmName = name,
				HasConverged = false,
				NumIterationsRequired = maxIterations,
				ResidualNormRatioEstimation = residualNormRatio
			};
		}

		private void CalculateAndPrintExactResidual()
		{
			var res = Vector.CreateZero(Rhs.Length);
			Matrix.Multiply(solution, res);
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
			/// Creates a new instance of <see cref="PcgAlgorithm"/>.
			/// </summary>
			public BlockPcgAlgorithm Build()
			{
				return new BlockPcgAlgorithm(BlockSize, ResidualTolerance, MaxIterationsProvider, Convergence, BlockResidualUpdater, 
					BetaCalculation);
			}
		}
	}
}
