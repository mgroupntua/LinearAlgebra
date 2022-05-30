using System;
using System.Collections.Generic;
using System.Text;
using MGroup.LinearAlgebra.Iterative;
using MGroup.LinearAlgebra.Iterative.PreconditionedConjugateGradient;
using MGroup.LinearAlgebra.Iterative.Termination;
using MGroup.LinearAlgebra.Distributed.IterativeMethods.Preconditioning;
using MGroup.MSolve.Solution.LinearSystem;

//TODOMPI: common IIterativeMethod interface for PCG, MINRES, GMRES. It is necessary so that the user of a DDM can choose the 
//  correct algorithm for his problem. 
//TODOMPI: generalize it to use IVector or even better an interface that only specifies the methods needed in iterative methods.
//TODOMPI: make it generic on the type of vector, which will be bounded by IVector or the dedicated interface (see previous TODO)
namespace MGroup.LinearAlgebra.Distributed.IterativeMethods.PCG
{
	public abstract class PcgAlgorithmBase : IDistributedIterativeMethod
	{
		protected readonly IPcgResidualConvergence convergence;
		protected readonly IPcgResidualUpdater residualUpdater;
		protected readonly bool throwIfNotConvergence;

		protected IGlobalVector direction;
		protected int iteration;
		protected IGlobalVector matrixTimesDirection;
		protected double paramBeta;
		protected IGlobalVector precondResidual;
		protected double resDotPrecondRes;
		protected double resDotPrecondResOld;
		protected IGlobalVector residual;
		protected IGlobalVector solution;
		protected double stepSize;

		protected PcgAlgorithmBase(double residualTolerance, IMaxIterationsProvider maxIterationsProvider,
			IPcgResidualConvergence convergence, IPcgResidualUpdater residualUpdater, bool throwIfNotConvergence)
		{
			this.ResidualTolerance = residualTolerance;
			this.MaxIterationsProvider = maxIterationsProvider;
			this.convergence = convergence;
			this.residualUpdater = residualUpdater;
			this.throwIfNotConvergence = throwIfNotConvergence;
		}

		public IMaxIterationsProvider MaxIterationsProvider { get; set; }

		public double ResidualTolerance { get; set; }

		/// <summary>
		/// The direction vector d, used to update the solution vector: x = x + α * d
		/// </summary>
		public IGlobalVector Direction => direction;

		/// <summary>
		/// The current iteration of the algorithm. It belongs to the interval [0, maxIterations).
		/// </summary>
		public int Iteration => iteration;

		/// <summary>
		/// The matrix A of the linear system or another object that implements matrix-vector multiplications.
		/// </summary>
		public MSolve.Solution.LinearSystem.ILinearTransformation Matrix { get; protected set; }

		/// <summary>
		/// The vector that results from <see cref="Matrix"/> * <see cref="Direction"/>.
		/// </summary>
		public IGlobalVector MatrixTimesDirection => matrixTimesDirection;

		/// <summary>
		/// The β parameter of Conjugate Gradient that ensures conjugacy between the direction vectors.
		/// </summary>
		public double ParamBeta => paramBeta;

		/// <summary>
		/// The preconditioner M, such that inv(M) ~= inv(A).
		/// </summary>
		public IPreconditioner Preconditioner { get; protected set; }

		/// <summary>
		/// The vector s = inv(M) * r
		/// </summary>
		public IGlobalVector PrecondResidual => precondResidual;

		/// <summary>
		/// The dot product r(t) * (inv(M) * r(t)) of the current iteration t.
		/// </summary>
		public double ResDotPrecondRes => resDotPrecondRes;

		/// <summary>
		/// The dot product r(t-1) * (inv(M) * r(t-1)) of the previous iteration t-1.
		/// </summary>
		public double ResDotPrecondResOld => resDotPrecondResOld;

		/// <summary>
		/// The residual vector r = b - A * x.
		/// </summary>
		public IGlobalVector Residual => residual;

		/// <summary>
		/// The right hand side of the linear system b = A * x.
		/// </summary>
		public IGlobalVector Rhs { get; protected set; }

		/// <summary>
		/// The current approximation to the solution of the linear system A * x = b
		/// </summary>
		public IGlobalVector Solution => solution;

		/// <summary>
		/// The step α taken along <see cref="Direction"/> to update the solution vector: x = x + α * d
		/// </summary>
		public double StepSize => stepSize;

		/// <summary>
		/// Releases references to the vectors and matrices used by this object and sets scalars to their default values.
		/// </summary>
		public virtual void Clear()
		{
			Matrix = null;
			Rhs = null;
			solution = null;
			residual = null;
			direction = null;
			matrixTimesDirection = null;
			Preconditioner = null;
			precondResidual = null;
			resDotPrecondRes = 0.0;
			resDotPrecondResOld = 0.0;
			stepSize = 0.0;
			paramBeta = 0.0;
			iteration = -1;
		}

		/// <summary>
		/// Solves the linear system A * x = b by solving the preconditioned system inv(P) * A * inv(P)^T * y = inv(P) * b, 
		/// where A = <paramref name="matrix"/>, b = <paramref name="rhsVector"/>, x is the solution, y = P^T * x,
		/// P*P^T = <paramref name="preconditioner"/>.
		/// Initially x = <paramref name="initialGuess"/> and then it converges to the solution.
		/// </summary>
		/// <param name="matrix">The matrix A of the linear system A * x = b. It must be symmetric positive definite.</param>
		/// <param name="rhs">
		/// The right hand side vector b of the linear system A * x = b. Constraints:
		/// <paramref name="rhs"/>.<see cref="IIndexable1D.Length"/> 
		/// == <paramref name="matrix"/>.<see cref="IIndexable2D.NumRows"/>.
		/// </param>
		/// <param name="preconditioner">
		/// A preconditioner matrix that is also symmetric positive definite and has the same dimensions as A.
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
		public virtual IterativeStatistics Solve(MSolve.Solution.LinearSystem.ILinearTransformation matrix, IPreconditioner preconditioner,
			IGlobalVector rhs, IGlobalVector solution, bool initialGuessIsZero) //TODO: find a better way to handle the case x0=0
		{
			this.Matrix = matrix;
			this.Preconditioner = preconditioner;
			this.Rhs = rhs;
			this.solution = solution;

			//TODOMPI: With distributed vectors/matrices, the dimensions may not be straightforward to calculate. In fact they
			//      may not be necessary in order to get the max iterations. E.g. In FETI methods, max iterations do not depend
			//      on the size of the global matrix of the interface problem, but are user defined usually.
			//int maxIterations = maxIterationsProvider.GetMaxIterations(matrix.NumColumns); 
			int maxIterations = ((FixedMaxIterationsProvider)MaxIterationsProvider).GetMaxIterations(-1);

			// r = b - A * x
			if (initialGuessIsZero) residual = rhs.Copy();
			else residual = ExactResidual.Calculate(matrix, rhs, solution);
			return SolveInternal(maxIterations, solution.CreateZero);

			//return Solve(new ExplicitMatrixTransformation(matrix), preconditioner, rhs, solution, initialGuessIsZero,
			//    zeroVectorInitializer);
		}

		protected abstract IterativeStatistics SolveInternal(int maxIterations, 
			Func<IGlobalVector> initializeZeroVector);
	}
}
