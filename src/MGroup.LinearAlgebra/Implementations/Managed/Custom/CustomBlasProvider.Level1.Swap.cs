namespace MGroup.LinearAlgebra.Implementations.Managed.Custom
{
	using System;
	using System.Diagnostics;

	public partial class CustomBlasProvider : IBlasProvider
	{
		public void Dswap(int n, double[] x, int offsetX, int incX, double[] y, int offsetY, int incY)
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
					double temp = x[baseX + ix];
					x[baseX + ix] = y[baseY + iy];
					y[baseY + iy] = temp;

					ix += incX;
					iy += incY;
				}

				return;
			}

			// Unit increments: loop unrolling (factor 3)
			int remainder = n % 3;

			// Cleanup loop
			for (int i = 1; i <= remainder; i++)
			{
				double temp = x[baseX + i];
				x[baseX + i] = y[baseY + i];
				y[baseY + i] = temp;
			}

			if (n < 3) return;

			// Main unrolled loop
			for (int i = remainder + 1; i <= n; i += 3)
			{
				double temp;

				temp = x[baseX + i];
				x[baseX + i] = y[baseY + i];
				y[baseY + i] = temp;

				temp = x[baseX + i + 1];
				x[baseX + i + 1] = y[baseY + i + 1];
				y[baseY + i + 1] = temp;

				temp = x[baseX + i + 2];
				x[baseX + i + 2] = y[baseY + i + 2];
				y[baseY + i + 2] = temp;
			}
		}
	}
}
