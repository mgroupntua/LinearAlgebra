namespace MGroup.LinearAlgebra.Implementations.Managed.Custom
{
	using System;
	using System.Diagnostics;

	public partial class CustomLapackProvider : ILapackProvider
	{
		// Applies a sequence of row interchanges (pivoting) to a matrix, typically as part of LU factorization.
		public void Dlaswp(int n, double[] a, int offsetA, int ldA, int k1, int k2, int[] ipiv, int offsetIpiv, int incx)
		{
			// Debug-only checks
			Debug.Assert(a != null);
			Debug.Assert(ipiv != null);
			Debug.Assert(ldA > 0);
			Debug.Assert(n >= 0);

			if (incx == 0) return;

			int ixStart, iStart, iEnd, iStep;

			if (incx > 0)
			{
				ixStart = k1;
				iStart = k1;
				iEnd = k2;
				iStep = 1;
			}
			else
			{
				ixStart = k1 + (k1 - k2) * incx;
				iStart = k2;
				iEnd = k1;
				iStep = -1;
			}

			// ------------------------------------------------------------
			// Process columns in blocks of 32 (loop unrolling optimization)
			// ------------------------------------------------------------
			int n32 = (n / 32) * 32;

			for (int j = 0; j < n32; j += 32)
			{
				int ix = ixStart;

				for (int i = iStart; i != iEnd + iStep; i += iStep)
				{
					int pivot = ipiv[offsetIpiv + ix - 1];

					if (pivot != i)
					{
						for (int jj = 0; jj < 32; jj++)
						{
							int col = j + jj;

							int index1 = offsetA + (i - 1) + col * ldA;
							int index2 = offsetA + (pivot - 1) + col * ldA;

							double temp = a[index1];
							a[index1] = a[index2];
							a[index2] = temp;
						}
					}

					ix += incx;
				}
			}

			// ------------------------------------------------------------
			// Remaining columns
			// ------------------------------------------------------------
			if (n32 < n)
			{
				int ix = ixStart;

				for (int i = iStart; i != iEnd + iStep; i += iStep)
				{
					int pivot = ipiv[offsetIpiv + ix - 1];

					if (pivot != i)
					{
						for (int col = n32; col < n; col++)
						{
							int index1 = offsetA + (i - 1) + col * ldA;
							int index2 = offsetA + (pivot - 1) + col * ldA;

							double temp = a[index1];
							a[index1] = a[index2];
							a[index2] = temp;
						}
					}

					ix += incx;
				}
			}
		}
	}
}
