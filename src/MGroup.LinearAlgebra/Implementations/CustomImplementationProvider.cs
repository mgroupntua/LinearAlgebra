namespace MGroup.LinearAlgebra.Implementations
{
	using System;
	using System.Collections.Generic;
	using System.Text;

	public class CustomImplementationProvider : IImplementationProvider
	{
		public CustomImplementationProvider(IBlasProvider blas, ISparseBlasProvider sparseBlas, ILapackProvider lapackProvider, 
			IReorderingProvider reordering) 
		{ 
			this.Blas = blas;
			this.SparseBlas = sparseBlas;
			this.LapackLinearEquations = new LapackLinearEquationsFacade(lapackProvider);
			this.LapackLeastSquares = new LapackLeastSquaresFacadeDouble(lapackProvider);
			this.LapackEigensystems = new LapackEigensystemsFacade(lapackProvider);
			this.Reordering = reordering;
		}

		public IBlasProvider Blas { get; }

		public LapackLinearEquationsFacade LapackLinearEquations { get; }

		public LapackLeastSquaresFacadeDouble LapackLeastSquares { get; }

		public LapackEigensystemsFacade LapackEigensystems { get; }

		public IReorderingProvider Reordering { get; }

		public ISparseBlasProvider SparseBlas { get; }
	}
}
