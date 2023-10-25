//using System;

//using MGroup.LinearAlgebra.Commons;
//using MGroup.LinearAlgebra.Exceptions;
//using MGroup.LinearAlgebra.Iterative.Termination;
//using MGroup.LinearAlgebra.Matrices;
//using MGroup.LinearAlgebra.Matrices.Operators;
//using MGroup.LinearAlgebra.Vectors;

//namespace MGroup.LinearAlgebra.Iterative.AlgebraicMultiGrid
//{

//	public class AlgebraicMultiGrid
//	{
//		private const string name = "Algebraic Multi-Grid";
//		private readonly double residualTolerance;
//		private readonly IMaxIterationsProvider maxIterationsProvider;

//		// <summary>
//		// The matrices used for prolongation operators
//		// </summary>
//		private readonly IMappingMatrix[] prolongationMatrixAtLevel;

//		// <summary>
//		// The matrices used for restriction operators
//		// </summary>
//		private readonly IMappingMatrix[] restrictionMatrixAtLevel;

//		private double resDotRes;
//		private double residualNormRatio;
//		private IVector solution;
//		private IVector residual;

//		// <summary>
//		// Internal iterative solver used for smoothing
//		// </summary>
//		private readonly GaussSeidel relaxationSolver;

//		// <summary>
//		// The offset position R_i:
//		// </summary>
//		private IVector[] offsetPositionVectors;

//		// <summary>
//		// The solution vectors x_i on level i, with x_0 being the main sought solution on level 0
//		// </summary>
//		private IVector[] solutionAtLevel;

//		// <summary>
//		// The rhs vectors b_i on level i, with b_0 being the main rhs on level 0
//		// </summary>
//		private IVector[] rhsAtLevel;

//		// <summary>
//		// The system matrices A_i for each level i, with A_0 being the main system matrix at level 0
//		// </summary>
//		private ILinearTransformation[] matrixAtLevel;

//		// <summary>
//		// The current convergence rate defined as solutionError^k/solutionError^{k-1}.
//		// The V-cycle iteration converges if convergenceRate_k < 1
//		// </summary>
//		private double convergenceRate;

//		// Properties //
//		/// <summary>
//		/// The matrix A of the linear system or another object that implements matrix-vector multiplications.
//		/// </summary>
//		public ILinearTransformation Matrix { get; private set; }

//		/// <summary>
//		/// The right hand side of the linear system b = A * x.
//		/// </summary>
//		public IVectorView Rhs { get; private set; }

//		/// <summary>
//		/// The current approximation to the solution of the linear system A * x = b
//		/// </summary>
//		public IVectorView Solution => solution;

//		/// <summary>
//		/// The current iteration of the algorithm. It belongs to the interval [0, maxIterations).
//		/// </summary>
//		public int Iteration { get; private set; }

//		/// <summary>
//		/// Implements the Algebraic Multi-Grid algorithm for solving linear systems.
//		/// Convergence is guaranteed only for strictly diagonally dominant or positive definite (symmetric) matrices.
//		/// Might converge in general matrix systems, as well, but with no guarantee.
//		/// Authors: Constantinos Atzarakis
//		/// </summary>
//		public AlgebraicMultiGrid(double residualTolerance, IMaxIterationsProvider maxIterationsProvider,
//								  GaussSeidel relaxationSolver, IMappingMatrix[] prolongationMatrixAtLevel,
//								  IMappingMatrix[] restrictionMatrixAtLevel)
//		{
//			if (residualTolerance <= 0)
//			{
//				throw new ArgumentOutOfRangeException(nameof(residualTolerance));
//			}

//			this.residualTolerance = residualTolerance;
//			this.maxIterationsProvider = maxIterationsProvider;
//			this.relaxationSolver = relaxationSolver;
//			this.prolongationMatrixAtLevel = prolongationMatrixAtLevel;
//			this.restrictionMatrixAtLevel = restrictionMatrixAtLevel;
//		}

//		public void Clear()
//		{
//			solutionAtLevel = null;
//			residual = null;
//			rhsAtLevel = null;
//			matrixAtLevel = null;
//			offsetPositionVectors = null;
//			Iteration = -1;
//			resDotRes = .0;
//		}

//		/// <summary>
//		/// Solves the linear system A * x = b, where A = <paramref name="matrix"/> and b = <paramref name="rhs"/>.
//		/// Initially x = <paramref name="initialGuess"/> and then it converges to the solution.
//		/// </summary>
//		/// <param name="matrix">The matrix A of the linear system A * x = b. It must be symmetric positive definite
//		/// or strictly diagonally dominant for ensured convergence.</param>
//		/// <param name="rhs">
//		/// The right hand side vector b of the linear system A * x = b. Constraints:
//		/// <paramref name="rhs"/>.<see cref="IIndexable1D.Length"/> 
//		/// == <paramref name="matrix"/>.<see cref="IIndexable2D.NumRows"/>.
//		/// </param>
//		/// <param name="solution">
//		/// The vector from which to start refining the solution vector x. Constraints:
//		/// <paramref name="solution"/>.<see cref="IIndexable1D.Length"/>
//		/// == <paramref name="matrix"/>.<see cref="IIndexable2D.NumColumns"/>.
//		/// </param>
//		/// <exception cref="NonMatchingDimensionsException">
//		/// Thrown if <paramref name="rhs"/> or <paramref name="solution"/> violate the described constraints.
//		/// </exception>
//		public IterativeStatistics Solve(ILinearTransformation matrix, IVectorView rhs, IVector solution)
//		{
//			//TODO: these will also be checked by the matrix vector multiplication.
//			Preconditions.CheckMultiplicationDimensions(matrix.NumColumns, solution.Length);
//			Preconditions.CheckSystemSolutionDimensions(matrix.NumRows, rhs.Length);

//			this.Matrix = matrix;
//			this.Rhs = rhs;
//			this.solution = solution;

//			residual = ExactResidual.Calculate(matrix, rhs, solution);

//			return SolveInternal(maxIterationsProvider.GetMaxIterations(matrix.NumColumns));
//		}

//		public IterativeStatistics Solve(IMatrixView matrix, IVectorView rhs, IVector solution)
//			=> Solve(new ExplicitMatrixTransformation(matrix), rhs, solution);

//		private IterativeStatistics SolveInternal(int maxiterations)
//		{
//			for (Iteration = 0; Iteration < maxiterations; Iteration++)
//			{
//				VCycleIteration();
//				residual = ExactResidual.Calculate(Matrix, Rhs, solutionAtLevel[0]);
//				resDotRes = residual.DotProduct(residual);
//				if (resDotRes < residualTolerance)
//				{
//					solution = solutionAtLevel[0];
//					return new IterativeStatistics
//						   {
//							   AlgorithmName = name, HasConverged = true, NumIterationsRequired = Iteration + 1,
//						   };
//				}
//			}

//			solution = solutionAtLevel[0];
//			return new IterativeStatistics
//				   {
//					   AlgorithmName = name, HasConverged = false, NumIterationsRequired = Iteration + 1,
//				   };
//		}

//		// <summary>
//		// Initializes the coarse levels of the multigrid hierarchy.
//		// Constructs A_i, b_i for each level i.
//		private void InitializeCoarseLevels()
//		{
//			throw new NotImplementedException();

//			matrixAtLevel[0] = Matrix;
//			rhsAtLevel[0] = Rhs as IVector;
//		}

//		// <summary>
//		// Performs a single V-cycle iteration of AMG.
//		// </summary>
//		private void VCycleIteration()
//		{
//			// Initialize solution and offset position vectors to 0.
//			// solutionAtLevel has one more element than offsetPositionVectors
//			solutionAtLevel[solutionAtLevel.Length].Clear();
//			for (int i = 0; i < solutionAtLevel.Length - 1; i++)
//			{
//				solutionAtLevel[i].Clear();
//				offsetPositionVectors[i].Clear();
//			}

//			for (int i = solutionAtLevel.Length; i >= 1; i--)
//			{
//				// Iterate on on level i, A_i * x_i = b_i
//				RelaxationAtLevel(i);

//				// apply prolongation operation
//				ProlongationMapping(i);
//			}

//			for (int i = 0; i < solutionAtLevel.Length - 1; i++)
//			{
//				// Iterate on on level i, A_i * x_i = b_i
//				RelaxationAtLevel(i);

//				// apply restriction operation
//				RestrictionMapping(i);
//			}
//		}


//		private void RelaxationAtLevel(int level)
//		{
//			// Iterate on on level i, A_i * x_i = b_i
//			relaxationSolver.Solve(matrixAtLevel[level], rhsAtLevel[level], solutionAtLevel[level]);
//		}

//		// <summary>
//		// Applys the restriction mapping and updates <see cref="offsetPositionVectors"/>
//		// and <see cref="rhsAtLevel"/> from level <paramref name="level"/> to level <paramref name="level"/> + 1.
//		// </summary>
//		private void RestrictionMapping(int level)
//		{
//			// Requires clustering info
//			// TODO: update x by interpolation matrix
//			// x_{i+1} = I_i^{i+1} x_i N^{i+1} ???
//			solutionAtLevel[level + 1] = restrictionMatrixAtLevel[level].Multiply(solutionAtLevel[level]);

//			// r_i = x_i - I^{i}_{i+1} x_{i+1}
//			offsetPositionVectors[level].AddIntoThis(solutionAtLevel[level]);
//			offsetPositionVectors[level].AxpyIntoThis(
//				restrictionMatrixAtLevel[level].Multiply(solutionAtLevel[level], offsetPositionVectors[level]),
//				-1.0);

//			// b_{i+1} = I^{i}_{i+1} (b_i - A_i r_i)
//			var diff = ExactResidual.Calculate(matrixAtLevel[level], rhsAtLevel[level], offsetPositionVectors[level]);
//			diff.ScaleIntoThis(-1.0);
//			rhsAtLevel[level + 1] = restrictionMatrixAtLevel[level].Multiply(diff);
//		}

//		private void ProlongationMapping(int level)
//		{
//			// apply prolongation operator X_i = P_i * X_{i+1} + R_i
//			prolongationMatrixAtLevel[level].Multiply(solutionAtLevel[level], solutionAtLevel[level - 1]);
//			solutionAtLevel[level - 1].AddIntoThis(offsetPositionVectors[level]);
//		}



//		/// <summary>
//		/// Constructs <see cref="AlgebraicMultiGrid"/> instances, allows the user to specify some or all of the required parameters and 
//		/// provides defaults for the rest.
//		/// Author: Serafeim Bakalakos
//		/// </summary>
//		public class Builder
//		{
//			/// <summary>
//			/// Specifies how to calculate the maximum iterations that the AGM algorithm will run for.
//			/// </summary>
//			public IMaxIterationsProvider MaxIterationsProvider { get; set; } =
//				new PercentageMaxIterationsProvider(1.0);

//			public double ResidualTolerance { get; set; } = 1E-10;

//			public GaussSeidel SmoothingSolver { get; set; } =
//				new GaussSeidel(1e-1, new PercentageMaxIterationsProvider(1.0));

//			/// <summary>
//			/// Creates a new instance of <see cref="AlgebraicMultiGrid"/>.
//			/// </summary>
//			public AlgebraicMultiGrid Build()
//				=> throw new NotImplementedException(); // AlgebraicMultiGrid(ResidualTolerance, MaxIterationsProvider, SmoothingSolver);
//		}
//	}
//}
