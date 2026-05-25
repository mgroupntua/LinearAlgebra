namespace MGroup.LinearAlgebra.Implementations.Managed.Custom
{
	using System;
	using System.Diagnostics;

	public partial class CustomBlasProvider : IBlasProvider
	{
		public void Daxpy(int n, double alpha, double[] x, int offsetX, int incX, double[] y, int offsetY, int incY)
		{
			Debug.Assert(x != null);
			Debug.Assert(y != null);
			Debug.Assert(offsetX >= 0);
			Debug.Assert(offsetY >= 0);
			Debug.Assert(incX != 0);
			Debug.Assert(incY != 0);

			if (n <= 0 || alpha == 0.0)
			{
				return;
			}

			if (incX != 1 || incY != 1)
			{
				int ix = offsetX;
				int iy = offsetY;

				if (incX < 0)
				{
					ix = offsetX + (n - 1) * incX;
				}

				if (incY < 0)
				{
					iy = offsetY + (n - 1) * incY;
				}

				for (int i = 0; i < n; i++)
				{
					y[iy] += alpha * x[ix];
					ix += incX;
					iy += incY;
				}

				return;
			}

			int remainder = n % 4;

			for (int i = 0; i < remainder; i++)
			{
				y[offsetY + i] += alpha * x[offsetX + i];
			}

			if (n < 4)
			{
				return;
			}

			for (int i = remainder; i < n; i += 4)
			{
				y[offsetY + i] += alpha * x[offsetX + i];
				y[offsetY + i + 1] += alpha * x[offsetX + i + 1];
				y[offsetY + i + 2] += alpha * x[offsetX + i + 2];
				y[offsetY + i + 3] += alpha * x[offsetX + i + 3];
			}
		}
	}
}
