namespace MGroup.LinearAlgebra.Implementations.Managed
{
	using System;
	using System.Diagnostics;

	public partial class ManagedBlasProvider : IBlasProvider
	{
		// Performs a rank-1 update: A = A + α*x*y^T
		public void Dger(int m, int n, double alpha, double[] x, int offsetX, int incX, double[] y, int offsetY, int incY, double[] a, int offsetA, int ldA)
		{
			// -----------------------------------------------------------
			// Debug-only argument validation (replaces XERBLA)
			// -----------------------------------------------------------
			Debug.Assert(m >= 0, "DGER: m must be >= 0");
			Debug.Assert(n >= 0, "DGER: n must be >= 0");
			Debug.Assert(incX != 0, "DGER: incX must not be 0");
			Debug.Assert(incY != 0, "DGER: incY must not be 0");
			Debug.Assert(ldA >= Math.Max(1, m), "DGER: ldA is invalid");
			Debug.Assert(x != null && y != null && a != null, "DGER: null array");

			// Quick return (BLAS behavior)
			if (m == 0 || n == 0 || alpha == 0.0)
				return;

			int yIndex = (incY > 0) ? 0 : (1 - n) * incY;

			if (incX == 1)
			{
				// Fast path: contiguous x
				for (int j = 0; j < n; j++)
				{
					double yj = y[offsetY + yIndex];

					if (yj != 0.0)
					{
						double scaledY = alpha * yj;
						int colOffset = offsetA + j * ldA;

						for (int i = 0; i < m; i++)
						{
							a[colOffset + i] += x[offsetX + i] * scaledY;
						}
					}

					yIndex += incY;
				}
			}
			else
			{
				int xStart = (incX > 0) ? 0 : (1 - m) * incX;

				for (int j = 0; j < n; j++)
				{
					double yj = y[offsetY + yIndex];

					if (yj != 0.0)
					{
						double scaledY = alpha * yj;

						int xIndex = xStart;
						int colOffset = offsetA + j * ldA;

						for (int i = 0; i < m; i++)
						{
							a[colOffset + i] += x[offsetX + xIndex] * scaledY;
							xIndex += incX;
						}
					}

					yIndex += incY;
				}
			}
		}
	}
}
