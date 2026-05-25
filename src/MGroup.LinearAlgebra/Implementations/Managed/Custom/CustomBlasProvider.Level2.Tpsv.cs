namespace MGroup.LinearAlgebra.Implementations.Managed.Custom
{
	using System;
	using System.Diagnostics;

	public partial class CustomBlasProvider : IBlasProvider
	{
		public void Dtpsv(StoredTriangle uplo, TransposeMatrix transA, DiagonalValues diag, int n,
			double[] a, int offsetA, double[] x, int offsetX, int incX)
		{
			bool unit = (diag == DiagonalValues.Unit) ? true : false;
			if (UseUpperImplementation(uplo, transA))
			{
				BackSubstitutionPackedColMajor(unit, n, a, offsetA, x, offsetX, incX);
			}
			else
			{
				ForwardSubstitutionPackedRowMajor(unit, n, a, offsetA, x, offsetX, incX);
			}
		}

		/// <summary>
		/// x = inv(U) * x
		/// </summary>
		internal static void BackSubstitutionPackedColMajor(bool unit, int n, double[] a, int offsetA,
			double[] x, int offsetX, int incX)
		{
			for (int i = n - 1; i >= 0; --i)
			{
				double dot = 0;
				for (int j = i + 1; j < n; ++j)
				{
					// index = i + ((j + 1) * j) / 2. They are not contigous during the dot products.
					dot += a[offsetA + i + ((j + 1) * j) / 2] * x[offsetX + j * incX];
				}

				// b[i] will not be used any more and can be replaced with x[i]
				if (unit) x[offsetX + i * incX] -= dot;
				else
				{
					int idxX = offsetX + i * incX;
					x[idxX] = (x[idxX] - dot) / a[offsetA + i + ((i + 1) * i) / 2];
				}
			}
		}

		/// <summary>
		/// x = inv(L) * x
		/// </summary>
		internal static void ForwardSubstitutionPackedRowMajor(bool unit, int n, double[] a, int offsetA,
			double[] x, int offsetX, int incX)
		{
			for (int i = 0; i < n; ++i)
			{
				double dot = 0;
				int rowOffset = offsetA + ((i + 1) * i) / 2;
				for (int j = 0; j < i; ++j)
				{
					// index = j + ((i + 1) * i) / 2. They are contigous during the dot products.
					dot += a[rowOffset + j] * x[offsetX + j * incX];
				}

				// b[i] will not be used any more and can be replaced with x[i]
				if (unit) x[offsetX + i * incX] -= dot;
				else
				{
					int idxX = offsetX + i * incX;
					x[idxX] = (x[idxX] - dot) / a[rowOffset + i];
				}
			}
		}
	}
}
