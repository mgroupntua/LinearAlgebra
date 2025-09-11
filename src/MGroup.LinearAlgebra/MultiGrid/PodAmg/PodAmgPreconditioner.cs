namespace MGroup.LinearAlgebra.AlgebraicMultiGrid.PodAmg
{
	using System;

	using MGroup.LinearAlgebra.AlgebraicMultiGrid;
	using MGroup.LinearAlgebra.Commons;
	using MGroup.LinearAlgebra.Exceptions;
	using MGroup.LinearAlgebra.Iterative.Preconditioning;
	using MGroup.LinearAlgebra.Matrices;
	using MGroup.LinearAlgebra.Triangulation;
	using MGroup.LinearAlgebra.Vectors;

	/// <summary>
	/// Preconditioner for iterative methods for linear systems, which uses the POD-AMG (Proper Orthogonal Decomposition - 
	/// Algebraic Multigrid) method.
	/// </summary>
	public class PodAmgPreconditioner : IPreconditioner
	{
		private readonly bool keepOnlyNonZeroPrincipalComponents;
		private readonly int numIterations;
		private readonly MultigridLevelSmoothing smoothing;

		private CsrMatrix fineMatrix;
		private CholeskyFull coarseMatrixFactorized;

		/// <summary>
		/// Prolongation/interpolation matrix. Restriction mastrix is assumed to be the transpose of this.
		/// </summary>
		private Matrix prolongation;

		/// <summary>
		/// Creates a new instance of <see cref="PodAmgAlgorithm"/> with the specified settings.
		/// </summary>
		/// <param name="keepOnlyNonZeroPrincipalComponents">How many principal components to keep during POD.</param>
		/// <param name="smoothing">
		/// Specifies the smoothing operator (e.g. Gauss-Seidel, Jacobi, SOR, ...) of the multigrid procedure.
		/// </param>
		/// <param name="numIterations">How many AMG cycles to perform each time the preconditioner is called.</param>
		public PodAmgPreconditioner(
			bool keepOnlyNonZeroPrincipalComponents, MultigridLevelSmoothing smoothing, int numIterations)
		{
			this.keepOnlyNonZeroPrincipalComponents = keepOnlyNonZeroPrincipalComponents;
			this.smoothing = smoothing;
			this.numIterations = numIterations;
		}

		/// <summary>
		/// Creates a new instance of <see cref="PodAmgPreconditioner"/> with the same settings as this instance.
		/// </summary>
		/// <returns>A new instance of <see cref="PodAmgPreconditioner"/>.</returns>
		public IPreconditioner CopyWithInitialSettings()
			=> new PodAmgPreconditioner(keepOnlyNonZeroPrincipalComponents, smoothing.CopyWithInitialSettings(), numIterations);

		/// <summary>
		/// Prepares the multigrid operators for the provided training data. This can be done only once, regardless of changes
		/// in the linear system matrix.
		/// </summary>
		/// <param name="sampleVectors">Matrix whose columns are the vectors to be used in POD.</param>
		/// <param name="numPrincipalComponents">How many principal components to keep during POD.</param>
		public void Initialize(Matrix sampleVectors, int numPrincipalComponents)
		{
			var pod = new ProperOrthogonalDecomposition(keepOnlyNonZeroPrincipalComponents);
			prolongation = pod.CalculatePrincipalComponents(sampleVectors.NumColumns, sampleVectors, numPrincipalComponents);
		}

		/// <summary>
		/// Solves the linear system of the preconditioning step.
		/// </summary>
		/// <param name="rhsVector">
		/// The right-hand-side vector of the preconditiong step. Usually the residual vector of the original linear system
		/// solver.
		/// </param>
		/// <param name="lhsVector">Initial guess for the solution vector. Usually it is zero.</param>
		public void SolveLinearSystem(IVectorView rhsVector, IVector lhsVector)
		{
			var rhs = (Vector)rhsVector;
			var solution = (Vector)lhsVector;
			solution.Clear();

			Preconditions.CheckSquareLinearSystemDimensions(fineMatrix, rhs, solution);
			var n0 = fineMatrix.NumRows;
			var r0 = Vector.CreateZero(n0);
			var e0 = Vector.CreateZero(n0);
			var n1 = coarseMatrixFactorized.Order;
			var r1 = Vector.CreateZero(n1);
			var e1 = Vector.CreateZero(n1);

			for (var i = 0; i < numIterations; i++)
			{
				// Pre-smoothing on lvl 0 to get an estimate of the solution x0. Use the x0 from previous cycles as initial guess.
				smoothing.ApplyPreSmoothers(rhs, solution);

				// Find the residual on lvl 0: r0=b-A0*x0
				//TODO: Use ExactResidual class for this
				fineMatrix.MultiplyIntoResult(solution, r0);
				r0.LinearCombinationIntoThis(-1.0, rhs, 1.0);

				// Restrict lvl 0 residual to lvl 1: r1 = P^T * r0
				prolongation.MultiplyIntoResult(r0, r1, transposeThis: true);

				// Find an estimate of the error on lvl 1 by solving exactly the system: A1*e1=r1.
				coarseMatrixFactorized.SolveLinearSystem(r1, e1);

				// Interpolate the lvl 1 error estimate to lvl 0: e0 = P * e1
				prolongation.MultiplyIntoResult(e1, e0, transposeThis: false);

				// Correct the solution estimate on lvl 0 using the interpolated error: x0 = x0 + e0
				solution.AddIntoThis(e0);

				// Post-smoothing on lvl 0 to further improve the solution estimate x0. Use the corrected x0 as initial guess.
				smoothing.ApplyPostSmoothers(rhs, solution);
			}
		}

		/// <summary>
		/// Prepares the multigrid operators for the provided matrix of the linear system. <see cref="Initialize(Matrix, int)"/> 
		/// must be called prior to this.
		/// </summary>
		/// <param name="matrix">The matrix of the original linear system in CSR format.</param>
		/// <param name="isPatternModified">
		/// True if the sparsity pattern of the matrix has been changed since the previous call of this method.
		/// </param>
		/// <exception cref="InvalidOperationException">
		/// Thrown if <see cref="Initialize(Matrix, int)"/> has not already been called.
		/// </exception>
		/// <exception cref="InvalidSparsityPatternException">
		/// Thrown if the matrix of the linear system is not in CSR format.
		/// </exception>
		public void UpdateMatrix(IMatrixView matrix, bool isPatternModified)
		{
			if (prolongation == null)
			{
				throw new InvalidOperationException("This preconditioner must be initialized first");
			}

			if (matrix is CsrMatrix csrMatrix)
			{
				fineMatrix = csrMatrix;
				smoothing.UpdateMatrix(fineMatrix, true);
				Matrix temp = fineMatrix.MultiplyRight(prolongation);
				Matrix coarseMatrix = prolongation.MultiplyRight(temp, transposeThis: true, transposeOther: false);
				coarseMatrixFactorized = CholeskyFull.Factorize(coarseMatrix.NumRows, coarseMatrix.RawData);
			}
			else
			{
				throw new InvalidSparsityPatternException("This preconditioner can be used only for matrices in CSR format");
			}
		}
	}
}
