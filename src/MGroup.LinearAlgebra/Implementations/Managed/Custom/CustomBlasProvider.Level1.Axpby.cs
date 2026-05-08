namespace MGroup.LinearAlgebra.Implementations.Managed.Custom
{
	using System;
	using System.Diagnostics;

	public partial class CustomBlasProvider : IBlasProvider
	{
		public void Daxpby(int n, double alpha, double[] x, int offsetX, int incX, double beta, double[] y, int offsetY, int incY)
		{
			Dscal(n, beta, y, offsetY, incY);
			Daxpy(n, alpha, x, offsetX, incX, y, offsetY, incY);
		}
	}
}
