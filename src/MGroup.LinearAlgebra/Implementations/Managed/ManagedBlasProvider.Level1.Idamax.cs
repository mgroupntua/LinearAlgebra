namespace MGroup.LinearAlgebra.Implementations.Managed
{
	using System;
	using System.Diagnostics;

	public partial class ManagedBlasProvider : IBlasProvider
	{
		public int Idamax(int n, double[] x, int offset, int incx)
		{
			if (n < 1 || incx <= 0)
			{
				return 0;
			}

			int maxIndex = 0;
			double maxValue = Math.Abs(x[offset]);

			if (incx == 1)
			{
				for (int i = 1; i < n; i++)
				{
					double value = Math.Abs(x[offset + i]);
					if (value > maxValue)
					{
						maxValue = value;
						maxIndex = i;
					}
				}
			}
			else
			{
				int idx = offset + incx;

				for (int i = 1; i < n; i++)
				{
					double value = Math.Abs(x[idx]);
					if (value > maxValue)
					{
						maxValue = value;
						maxIndex = i;
					}

					idx += incx;
				}
			}

			return maxIndex;
		}
	}
}
