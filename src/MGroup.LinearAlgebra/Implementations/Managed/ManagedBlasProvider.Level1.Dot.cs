namespace MGroup.LinearAlgebra.Implementations.Managed
{
	using System;
	using System.Diagnostics;

	public partial class ManagedBlasProvider : IBlasProvider
	{
		public double Ddot(int n, double[] x, int offsetX, int incX, double[] y, int offsetY, int incY)
		{
			Debug.Assert(x != null);
			Debug.Assert(y != null);
			Debug.Assert(offsetX >= 0);
			Debug.Assert(offsetY >= 0);
			Debug.Assert(incX != 0);
			Debug.Assert(incY != 0);

			double result = 0.0;

			if (n <= 0)
			{
				return result;
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
					result += x[ix] * y[iy];
					ix += incX;
					iy += incY;
				}

				return result;
			}

			int remainder = n % 5;
			for (int i = 0; i < remainder; i++)
			{
				result += x[offsetX + i] * y[offsetY + i];
			}

			if (n < 5)
			{
				return result;
			}

			for (int i = remainder; i < n; i += 5)
			{
				result += x[offsetX + i] * y[offsetY + i]
						+ x[offsetX + i + 1] * y[offsetY + i + 1]
						+ x[offsetX + i + 2] * y[offsetY + i + 2]
						+ x[offsetX + i + 3] * y[offsetY + i + 3]
						+ x[offsetX + i + 4] * y[offsetY + i + 4];
			}

			return result;
		}
	}
}
