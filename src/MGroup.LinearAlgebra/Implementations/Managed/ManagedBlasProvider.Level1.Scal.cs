namespace MGroup.LinearAlgebra.Implementations.Managed
{
	using System;
	using System.Diagnostics;

	public partial class ManagedBlasProvider : IBlasProvider
	{
		public void Dscal(int n, double alpha, double[] x, int offsetX, int incX)
		{
			Debug.Assert(x != null);
			Debug.Assert(offsetX >= 0);
			Debug.Assert(incX != 0);

			if (n <= 0 || incX <= 0)
			{
				return;
			}

			int start = offsetX;

			// Case 1: General increment (incX != 1)
			if (incX != 1)
			{
				int index = start;
				for (int i = 0; i < n; i++)
				{
					x[index] *= alpha;
					index += incX;
				}

				return;
			}

			// Case 2: Unit increment (incX == 1)
			// Uses loop unrolling for better performance
			int remainder = n % 5;

			// Cleanup loop (process elements not divisible by 5)
			for (int i = 0; i < remainder; i++)
			{
				x[start + i] *= alpha;
			}

			if (n < 5)
			{
				return;
			}

			// Main unrolled loop (process 5 elements per iteration)
			for (int i = remainder; i < n; i += 5)
			{
				x[start + i] *= alpha;
				x[start + i + 1] *= alpha;
				x[start + i + 2] *= alpha;
				x[start + i + 3] *= alpha;
				x[start + i + 4] *= alpha;
			}
		}
	}
}
