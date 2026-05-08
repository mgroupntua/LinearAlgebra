namespace MGroup.LinearAlgebra.Implementations.Managed.Custom
{
	using System;
	using System.Diagnostics;

	public partial class CustomLapackProvider : ILapackProvider
	{
		// DGETF2 performs an unblocked LU factorization with partial pivoting of a general matrix.
		public void Dgetf2(int m, int n, double[] a, int offsetA, int ldA, int[] ipiv, int offsetIpiv, out int info)
		{
			info = 0;

			// ---------------------------
			// Argument checks
			// ---------------------------
			if (m < 0) info = -1;
			else if (n < 0) info = -2;
			else if (ldA < Math.Max(1, m)) info = -4;

			if (info != 0)
			{
				string paramName = info switch
				{
					-1 => nameof(m),
					-2 => nameof(n),
					-4 => nameof(ldA),
					_ => "unknown"
				};

				throw new ArgumentException($"DGETF2: invalid value for parameter '{paramName}'.");
			}

			if (m == 0 || n == 0)
				return;

			// Machine safe minimum (used to avoid overflow)
			double safeMin = Dlamch(DlamchParam.SafeMinimum);

			int minMN = Math.Min(m, n);

			// ---------------------------
			// Main loop over columns
			// ---------------------------
			for (int j = 0; j < minMN; j++)
			{
				// ---------------------------
				// Pivot selection
				// ---------------------------
				int pivotOffset = offsetA + j + j * ldA;
				int pivotIndex = j + blas.Idamax(m - j, a, pivotOffset, 1);

				ipiv[offsetIpiv + j] = pivotIndex;

				double pivotValue = a[offsetA + pivotIndex + j * ldA];

				if (pivotValue != 0.0)
				{
					// ---------------------------
					// Row interchange
					// ---------------------------
					if (pivotIndex != j)
					{
						blas.Dswap(n,
							  a, offsetA + j, ldA,
							  a, offsetA + pivotIndex, ldA);
					}

					// ---------------------------
					// Compute multipliers
					// ---------------------------
					if (j < m - 1)
					{
						double diag = a[offsetA + j + j * ldA];

						if (Math.Abs(diag) >= safeMin)
						{
							// Use BLAS scaling
							blas.Dscal(m - j - 1,
								  1.0 / diag,
								  a, offsetA + (j + 1) + j * ldA, 1);
						}
						else
						{
							// Safe manual division (avoids overflow)
							for (int i = j + 1; i < m; i++)
							{
								a[offsetA + i + j * ldA] /= diag;
							}
						}
					}
				}
				else if (info == 0)
				{
					// Singular matrix
					info = j + 1;
				}

				// ---------------------------
				// Trailing submatrix update
				// ---------------------------
				if (j < minMN - 1)
				{
					blas.Dger(m - j - 1, n - j - 1, -1.0,
						a, offsetA + (j + 1) + j * ldA, 1,
						a, offsetA + j + (j + 1) * ldA, ldA,
						a, offsetA + (j + 1) + (j + 1) * ldA, ldA);
				}
			}
		}
	}
}
