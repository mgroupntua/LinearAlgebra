
namespace MGroup.LinearAlgebra.Implementations.NativeWin64
{
	using System;
	using System.Collections.Generic;
	using System.Text;

	using MGroup.LinearAlgebra.Implementations.NativeWin64.MKL;

	/// <summary>
	/// Providers that call the highly optimized Intel Math Kernel Library (native dlls). Note that Intel MKL relies heavily
	/// on OpenMP (and sometimes TBB), which might not be desirable if the user also wants to fully utilize the available
	/// CPU cores in their own multi-threaded code.
	/// </summary>
	public class NativeWin64ImplementationProvider : IImplementationProvider
	{
		public NativeWin64ImplementationProvider()
		{
			Blas = MklBlasProvider.UniqueInstance;
			SparseBlas = MklSparseBlasProvider.UniqueInstance;
			LapackLinearEquations = new LapackLinearEquationsFacade(MklLapackProvider.UniqueInstance);
			LapackLeastSquares = new LapackLeastSquaresFacadeDouble(MklLapackProvider.UniqueInstance);
			LapackEigensystems = new LapackEigensystemsFacade(MklLapackProvider.UniqueInstance);
			Reordering = new SuiteSparseReorderingProvider();
		}

		public IBlasProvider Blas { get; }

		public LapackLinearEquationsFacade LapackLinearEquations { get; }

		public LapackLeastSquaresFacadeDouble LapackLeastSquares { get; }

		public LapackEigensystemsFacade LapackEigensystems { get; }

		public IReorderingProvider Reordering { get; }

		public ISparseBlasProvider SparseBlas { get; }
	}
}
