namespace MGroup.LinearAlgebra.Implementations.Managed
{
	using System;
	using System.Diagnostics;

	public partial class ManagedBlasProvider : IBlasProvider
	{
		public void Dspmv(StoredTriangle uplo, int n,
			double alpha, double[] a, int offsetA, double[] x, int offsetX, int incX,
			double beta, double[] y, int offsetY, int incY)
		{
			// y = alpha * L * x + beta * y 
			LowerTimesVectorPackedRowMajor(
				Diagonal.Regular, n, alpha, a, offsetA, x, offsetX, incX, beta, y, offsetY, incY);

			// y = alpha * U * x + y, where U has 0 diagonal
			UpperTimesVectorPackedColMajor(
				Diagonal.Zero, n, alpha, a, offsetA, x, offsetX, incX, 1.0, y, offsetY, incY);
		}
	}
}
