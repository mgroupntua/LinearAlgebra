namespace MGroup.LinearAlgebra.Implementations
{
	using System;
	using System.Collections.Generic;
	using System.Text;

	public class CustomImplementationProvider : IImplementationProvider
	{
		public CustomImplementationProvider(IBlasProvider blas, ISparseBlasProvider sparseBlas, ILapackProvider lapackProvider) 
		{ 
			this.Blas = blas;
			this.SparseBlas = sparseBlas;
			this.LapackLinearEquations = new LapackLinearEquationsFacade(lapackProvider);
			this.LapackLeastSquares = new LapackLeastSquaresFacadeDouble(lapackProvider);
			this.LapackEigensystems = new LapackEigensystemsFacade(lapackProvider);
		}

		public IBlasProvider Blas { get; }

		public LapackLinearEquationsFacade LapackLinearEquations { get; }

		public LapackLeastSquaresFacadeDouble LapackLeastSquares { get; }

		public LapackEigensystemsFacade LapackEigensystems { get; }

		public ISparseBlasProvider SparseBlas { get; }
	}
}
