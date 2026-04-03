namespace MGroup.LinearAlgebra.Implementations.Managed
{
	using System;
	using System.Diagnostics;

	using DotNumerics.LinearAlgebra.CSLapack;

	public partial class ManagedBlasProvider : IBlasProvider
	{
		public double Dasum(int n, double[] x, int offsetX, int incX)
		{
			if (n <= 0 || incX <= 0) return 0.0;

			int baseX = offsetX - 1;
			double sum = 0.0;

			// General increments
			if (incX != 1)
			{
				int ix = 1;

				for (int i = 0; i < n; i++)
				{
					sum += Math.Abs(x[baseX + ix]);
					ix += incX;
				}

				return sum;
			}

			// Unit increments: loop unrolling (factor 6)
			int remainder = n % 6;

			// Cleanup loop
			for (int i = 1; i <= remainder; i++)
			{
				sum += Math.Abs(x[baseX + i]);
			}

			if (n < 6) return sum;

			// Main unrolled loop
			for (int i = remainder + 1; i <= n; i += 6)
			{
				sum += Math.Abs(x[baseX + i])
					 + Math.Abs(x[baseX + i + 1])
					 + Math.Abs(x[baseX + i + 2])
					 + Math.Abs(x[baseX + i + 3])
					 + Math.Abs(x[baseX + i + 4])
					 + Math.Abs(x[baseX + i + 5]);
			}

			return sum;
		}
	}
}
