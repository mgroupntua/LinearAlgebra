namespace MGroup.LinearAlgebra.Implementations.Managed
{
	using System;
	using System.Diagnostics;

	public partial class ManagedBlasProvider : IBlasProvider
	{
		public void Dcopy(int n, double[] x, int offsetX, int incX, double[] y, int offsetY, int incY)
		{
			if (n <= 0) return;

			int baseX = offsetX - 1;
			int baseY = offsetY - 1;

			// General increments
			if (incX != 1 || incY != 1)
			{
				int ix = 1;
				int iy = 1;

				if (incX < 0) ix = (1 - n) * incX + 1;
				if (incY < 0) iy = (1 - n) * incY + 1;

				for (int i = 0; i < n; i++)
				{
					y[baseY + iy] = x[baseX + ix];
					ix += incX;
					iy += incY;
				}

				return;
			}

			// Unit increments: loop unrolling (factor 7)
			int remainder = n % 7;

			// Cleanup loop
			for (int i = 1; i <= remainder; i++)
			{
				y[baseY + i] = x[baseX + i];
			}

			if (n < 7) return;

			// Main unrolled loop
			for (int i = remainder + 1; i <= n; i += 7)
			{
				y[baseY + i] = x[baseX + i];
				y[baseY + i + 1] = x[baseX + i + 1];
				y[baseY + i + 2] = x[baseX + i + 2];
				y[baseY + i + 3] = x[baseX + i + 3];
				y[baseY + i + 4] = x[baseX + i + 4];
				y[baseY + i + 5] = x[baseX + i + 5];
				y[baseY + i + 6] = x[baseX + i + 6];
			}
		}
	}
}
