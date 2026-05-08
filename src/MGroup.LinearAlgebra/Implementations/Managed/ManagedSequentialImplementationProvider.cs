namespace MGroup.LinearAlgebra.Implementations.Managed
{
	using System;
	using System.Collections.Generic;
	using System.Text;

	using MGroup.LinearAlgebra.Implementations.Managed.Custom;
	using MGroup.LinearAlgebra.Triangulation;

	/// <summary>
	/// Providers that call one or more linear algebra libraries written in managed C# code. Those libraries could be 3rd 
	/// party, custom ones or any combination thereof. They should be used if the libraries required by the other providers 
	/// are not installed on the user's system.
	/// </summary>
	public class ManagedSequentialImplementationProvider : IImplementationProvider
	{
		private readonly double luPivotTolerance;

		public ManagedSequentialImplementationProvider(double luPivotTolerance = LUCSparseNet.DefaultPivotTolerance)
		{
			this.luPivotTolerance = luPivotTolerance;
			Blas = CustomBlasProvider.UniqueInstance;
			SparseBlas = ManagedSparseBlasProvider.UniqueInstance;
			LapackLinearEquations = new LapackLinearEquationsFacade(DotNumericsLapackProvider.UniqueInstance);
			LapackLeastSquares = new LapackLeastSquaresFacadeDouble(DotNumericsLapackProvider.UniqueInstance);
			LapackEigensystems = new LapackEigensystemsFacade(DotNumericsLapackProvider.UniqueInstance);
			Reordering = new ManagedReorderingProvider();
		}

		public IBlasProvider Blas { get; }

		public LapackLinearEquationsFacade LapackLinearEquations { get; }

		public LapackLeastSquaresFacadeDouble LapackLeastSquares { get; }

		public LapackEigensystemsFacade LapackEigensystems { get; }

		public IReorderingProvider Reordering { get; }

		public ISparseBlasProvider SparseBlas { get; }

		public ILUCscFactorization CreateLUTriangulation() => new LUCSparseNet(luPivotTolerance);

		public ICholeskySymmetricCsc CreateCholeskyTriangulation() => new CholeskyCSparseNet();
	}
}
