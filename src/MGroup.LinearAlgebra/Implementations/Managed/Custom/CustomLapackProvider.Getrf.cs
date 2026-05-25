namespace MGroup.LinearAlgebra.Implementations.Managed.Custom
{
	using System;
	using System.Diagnostics;

	//using DotNumerics.LinearAlgebra.CSLapack;

	public partial class CustomLapackProvider : ILapackProvider
	{
		/// <summary>
		/// Computes LU factorization of a general MxN matrix using partial pivoting.
		/// </summary>
		public void Dgetrf(int m, int n, double[] a, int offsetA, int ldA, int[] ipiv, int offsetIpiv, ref int info)
		{
			throw new Exception("The tests do not pass with the cleaned up code. Use CSLapack.DGETRF instead"); 
			//dgetrf.Run(m, n, ref a, offsetA, ldA, ref ipiv, offsetIpiv, ref info);
			//return;

			info = 0;

			// ---------------------------
			// Argument validation (debug only)
			// ---------------------------
			Debug.Assert(m >= 0, "M must be >= 0");
			Debug.Assert(n >= 0, "N must be >= 0");
			Debug.Assert(ldA >= Math.Max(1, m), "ldA is invalid");

			if (m == 0 || n == 0)
			{
				return;
			}

			int minMN = Math.Min(m, n);

			// ---------------------------
			// Determine block size
			// ---------------------------
			int nb = Ilaenv(IlaenvSpec.OptimalBlockSize, LapackDataType.DoubleReal, LapackMatrixType.General, LapackOperation.TRF, LapackOptions.None,
				m, n, -1, -1);

			// ---------------------------
			// Use unblocked code
			// ---------------------------
			if (nb <= 1 || nb >= minMN)
			{
				Dgetf2(m, n, a, offsetA, ldA, ipiv, offsetIpiv, out info);
				return;
			}

			// ---------------------------
			// Blocked algorithm
			// ---------------------------
			for (int j = 0; j < minMN; j += nb)
			{
				int jb = Math.Min(minMN - j, nb);

				int panelOffset = offsetA + j + j * ldA;

				// Factor current panel
				Dgetf2(m - j, jb, a, panelOffset, ldA, ipiv, offsetIpiv + j, out int iinfo);

				if (info == 0 && iinfo > 0)
				{
					info = iinfo + j;
				}

				// Adjust pivot indices
				for (int i = j; i < j + jb && i < m; i++)
				{
					ipiv[offsetIpiv + i] += j;
				}

				// Apply row interchanges to left block
				Dlaswp(j, a, offsetA, ldA, j, j + jb - 1, ipiv, offsetIpiv, 1);

				if (j + jb < n)
				{
					// Apply row interchanges to right block
					Dlaswp(n - j - jb, a, offsetA + (j + jb) * ldA, ldA, j, j + jb - 1, ipiv, offsetIpiv, 1);

					// Compute block row of U
					blas.Dtrsm(MultiplicationSide.Left, StoredTriangle.Lower, TransposeMatrix.NoTranspose, DiagonalValues.Unit, jb, n - j - jb, 1.0, a, panelOffset, ldA, a, offsetA + j + (j + jb) * ldA, ldA);

					if (j + jb < m)
					{
						// Update trailing submatrix: A22 = A22 - L21 * U12
						// where L21 = subdiagonal block (below panel),
						// U12 = block row to the right,
						// A22 = bottom - right trailing matrix
						int mTrailing = m - j - jb;
						int nTrailing = n - j - jb;
						int kBlock = jb;
						int offsetL21 = offsetA + (j + jb) + j * ldA;             // L21 block
						int offsetU12 = offsetA + j + (j + jb) * ldA;             // U12 block
						int offsetA22 = offsetA + (j + jb) + (j + jb) * ldA;      // A22 block

						blas.Dgemm(TransposeMatrix.NoTranspose, TransposeMatrix.NoTranspose, mTrailing, nTrailing, kBlock, -1.0, a, offsetL21, ldA, a, offsetU12, ldA, 1.0, a, offsetA22, ldA);
					}
				}
			}
		}
	}
}
