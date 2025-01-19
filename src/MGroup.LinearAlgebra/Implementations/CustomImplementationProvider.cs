namespace MGroup.LinearAlgebra.Implementations
{
	using System;
	using System.Collections.Generic;
	using System.Text;

	using MGroup.LinearAlgebra.Triangulation;

	public class CustomImplementationProvider : IImplementationProvider
	{
		private readonly Func<ICholeskySymmetricCsc> createCholeskyTriangulation;
		private readonly Func<ILUCscFactorization> createLUTriangulation;

		public CustomImplementationProvider(IBlasProvider blas, ISparseBlasProvider sparseBlas, ILapackProvider lapackProvider,
			IReorderingProvider reordering,
			Func<ILUCscFactorization> createLUTriangulation,
			Func<ICholeskySymmetricCsc> createCholeskyTriangulation)
		{ 
			this.Blas = blas;
			this.SparseBlas = sparseBlas;
			this.LapackLinearEquations = new LapackLinearEquationsFacade(lapackProvider);
			this.LapackLeastSquares = new LapackLeastSquaresFacadeDouble(lapackProvider);
			this.LapackEigensystems = new LapackEigensystemsFacade(lapackProvider);
			this.Reordering = reordering;
			this.createLUTriangulation = createLUTriangulation;
			this.createCholeskyTriangulation = createCholeskyTriangulation;
		}

		public IBlasProvider Blas { get; }

		public LapackLinearEquationsFacade LapackLinearEquations { get; }

		public LapackLeastSquaresFacadeDouble LapackLeastSquares { get; }

		public LapackEigensystemsFacade LapackEigensystems { get; }

		public IReorderingProvider Reordering { get; }

		public ISparseBlasProvider SparseBlas { get; }

		public ILUCscFactorization CreateLUTriangulation() => createLUTriangulation();

		public ICholeskySymmetricCsc CreateCholeskyTriangulation() => createCholeskyTriangulation();
	}
}
