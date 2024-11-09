namespace MGroup.LinearAlgebra.Implementations.Managed
{
	using System;
	using System.Collections.Generic;
	using System.Text;

	/// <summary>
	/// Providers that call one or more linear algebra libraries written in managed C# code. Those libraries could be 3rd 
	/// party, custom ones or any combination thereof. They should be used if the libraries required by the other providers 
	/// are not installed on the user's system.
	/// </summary>
	public class ManagedSequentialImplementationProvider : IImplementationProvider
	{
		public ManagedSequentialImplementationProvider()
		{
			Blas = ManagedBlasProvider.UniqueInstance;
			SparseBlas = ManagedSparseBlasProvider.UniqueInstance;
			LapackLinearEquations = new LapackLinearEquationsFacade(ManagedLapackProvider.UniqueInstance);
			LapackLeastSquares = new LapackLeastSquaresFacadeDouble(ManagedLapackProvider.UniqueInstance);
			LapackEigensystems = new LapackEigensystemsFacade(ManagedLapackProvider.UniqueInstance);
		}

		public IBlasProvider Blas { get; }

		public LapackLinearEquationsFacade LapackLinearEquations { get; }

		public LapackLeastSquaresFacadeDouble LapackLeastSquares { get; }

		public LapackEigensystemsFacade LapackEigensystems { get; }

		public ISparseBlasProvider SparseBlas { get; }
	}
}
