namespace MGroup.LinearAlgebra.Implementations
{
	using System;
	using System.Collections.Generic;
	using System.Text;

	public interface IImplementationProvider
	{
		public IBlasProvider Blas { get; }

		public LapackLinearEquationsFacade LapackLinearEquations { get; }

		public LapackLeastSquaresFacadeDouble LapackLeastSquares { get; }

		public LapackEigensystemsFacade LapackEigensystems { get; }

		public IReorderingProvider Reordering { get; }

		public ISparseBlasProvider SparseBlas { get; }
	}
}
