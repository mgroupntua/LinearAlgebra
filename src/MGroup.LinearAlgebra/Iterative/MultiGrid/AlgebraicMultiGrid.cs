using System;
using MGroup.LinearAlgebra.Iterative.Termination;
using MGroup.LinearAlgebra.Vectors;
using MGroup.LinearAlgebra.Matrices;


namespace MGroup.LinearAlgebra.Iterative.MultiGrid;

public class AlgebraicMultiGrid
{
	private const string name = "Algebraic Multi-Grid";
	private readonly double residualTolerance;
	private readonly IMaxIterationsProvider maxIterationsProvider;

	private IVector solution;
	private IVector residual;

	/// <summary>
	/// Implements the Algebraic Multi-Grid algorithm for solving linear systems.
	/// Convergence is guaranteed only for strictly diagonally dominant or positive definite (symmetric) matrices.
	/// Might converge in general matrix systems, as well, but with no guarantee.
	/// Authors: Constantinos Atzarakis
	/// </summary>
	public AlgebraicMultiGrid(double residualTolerance, IMaxIterationsProvider maxIterationsProvider)
	{
		if (residualTolerance <= 0)
		{
			throw new ArgumentOutOfRangeException(nameof(residualTolerance));
		}

		this.residualTolerance = residualTolerance;
		this.maxIterationsProvider = maxIterationsProvider;
	}

	/// <summary>
	/// The current iteration of the algorithm. It belongs to the interval [0, maxIterations).
	/// </summary>	
	public int Iteration { get; private set; }

	/// <summary>
	/// The matrix A of the linear system or another object that implements matrix-vector multiplications.
	/// </summary>
	public ILinearTransformation Matrix { get; private set; }

	/// <summary>
	/// The current approximation to the solution of the linear system A * x = b
	/// </summary>
	public IVectorView Solution => solution;

	/// <summary>
	/// Solves the linear system A * x = b, where A = <paramref name="matrix"/> and b = <paramref name="rhs"/>.
	/// Initially x = <paramref name="initialGuess"/> and then it converges to the solution.
	/// </summary>
	/// <param name="matrix">The matrix A of the linear system A * x = b. It must be symmetric positive definite
	/// or strictly diagonally dominant for ensured convergence.</param>
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
	public IterativeStatistics Solve(ILinearTransformation matrix, IVectorView rhs, IVectorView solution)
	{
		throw new NotImplementedException();
	}

	public IterativeStatistics Solve(IMatrixView matrix, IVectorView rhs, IVector solution)
		=> Solve(new ExplicitMatrixTransformation(matrix), rhs, solution);

	public void Clear()
	{
		throw new NotImplementedException();
	}

	/// <summary>
	/// Constructs <see cref="AlgebraicMultiGrid"/> instances, allows the user to specify some or all of the required parameters and 
	/// provides defaults for the rest.
	/// Author: Serafeim Bakalakos
	/// </summary>
	public class Builder
	{
		/// <summary>
		/// Specifies how to calculate the maximum iterations that the AGM algorithm will run for.
		/// </summary>
		public IMaxIterationsProvider MaxIterationsProvider { get; set; } = new PercentageMaxIterationsProvider(1.0);

		public double ResidualTolerance { get; set; } = 1E-10;

		/// <summary>
		/// Creates a new instance of <see cref="AlgebraicMultiGrid"/>.
		/// </summary>
		public AlgebraicMultiGrid Build()
			=> new AlgebraicMultiGrid(ResidualTolerance, MaxIterationsProvider);
	}

}
