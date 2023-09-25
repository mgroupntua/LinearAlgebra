namespace MGroup.LinearAlgebra.Iterative.GaussSeidel
{
	using System;
	using System.Collections.Generic;
	using System.Text;

	using MGroup.LinearAlgebra.Matrices;
	using MGroup.LinearAlgebra.Vectors;

	public interface IGaussSeidelIteration : IDisposable
	{
		void Initialize(IMatrixView matrix);

		/// <summary>
		/// Performs 1 backward Gauss-Seidel iteration.
		/// </summary>
		/// <param name="rhsVector">Will not check if the dimensions are compatible with the matrix.</param>
		/// <param name="lhsVector">Will not check if the dimensions are compatible with the matrix.</param>
		void GaussSeidelBackwardIteration(IVectorView rhsVector, IVector lhsVector);

		/// <summary>
		/// Performs 1 forward Gauss-Seidel iteration.
		/// </summary>
		/// <param name="rhsVector">Will not check if the dimensions are compatible with the matrix.</param>
		/// <param name="lhsVector">Will not check if the dimensions are compatible with the matrix.</param>
		void GaussSeidelForwardIteration(IVectorView rhsVector, IVector lhsVector);
	}

	public interface IGaussSeidelIterationBuilder
	{
		IGaussSeidelIteration Create();
	}
}
