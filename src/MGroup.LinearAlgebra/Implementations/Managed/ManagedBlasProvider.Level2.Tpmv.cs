namespace MGroup.LinearAlgebra.Implementations.Managed
{
	using System;
	using System.Diagnostics;

	public partial class ManagedBlasProvider : IBlasProvider
	{
		public void Dtpmv(StoredTriangle uplo, TransposeMatrix transA, DiagonalValues diag, int n, double[] a, int offsetA, double[] x, int offsetX, int incX)
		{
			// The copy may be avoidable in trangular operations, if we start the dot products from the bottom
			var input = new double[x.Length];
			Array.Copy(x, input, x.Length);

			Diagonal managedDiag = (diag == DiagonalValues.NonUnit) ?
				Diagonal.Regular : Diagonal.Unit;
			if (UseUpperImplementation(uplo, transA))
			{
				UpperTimesVectorPackedColMajor(
					managedDiag, n, 1.0, a, offsetA, input, offsetX, incX, 0.0, x, offsetX, incX);
			}
			else
			{
				LowerTimesVectorPackedRowMajor(
					managedDiag, n, 1.0, a, offsetA, input, offsetX, incX, 0.0, x, offsetX, incX);
			}
		}

		/// <summary>
		/// y = alpha * A * x + beta * y
		/// </summary>
		internal static void LowerTimesVectorPackedRowMajor(Diagonal diag, int n, double alpha, double[] a, int offsetA,
			double[] x, int offsetX, int incX, double beta, double[] y, int offsetY, int incY)
		{
			for (int i = 0; i < n; ++i)
			{
				// Row * vector without the diagonal
				double dot = 0.0;
				int rowOffset = offsetA + ((i + 1) * i) / 2;
				for (int j = 0; j < i; ++j)
				{
					// index = j + ((i + 1) * i) / 2. They are contigous during the dot products.
					dot += a[rowOffset + j] * x[offsetX + j * incX];
				}

				// Take into account the diagonal
				if (diag == Diagonal.Regular) dot += a[rowOffset + i] * x[offsetX + i * incX];
				else if (diag == Diagonal.Unit) dot += x[offsetX + i * incX];
				// else the diagonal is ignored altogether (we operate only on the superdiagonal part)

				// Take into account the rest and store it
				int idxY = offsetY + i * incY;
				y[idxY] = alpha * dot + beta * y[idxY];
			}
		}

		/// <summary>
		/// y = alpha * A * x + beta * y
		/// </summary>
		internal static void UpperTimesVectorPackedColMajor(Diagonal diag, int n, double alpha, double[] a, int offsetA,
			double[] x, int offsetX, int incX, double beta, double[] y, int offsetY, int incY)
		{
			for (int i = 0; i < n; ++i)
			{
				// Row * vector without the diagonal
				double dot = 0.0;
				for (int j = i + 1; j < n; ++j)
				{
					// index = i + ((j + 1) * j) / 2. They are not contigous during the dot products.
					dot += a[i + ((j + 1) * j) / 2] * x[offsetX + j * incX];
				}

				// Take into account the diagonal
				if (diag == Diagonal.Regular) dot += a[offsetA + i + ((i + 1) * i) / 2] * x[offsetX + i * incX];
				else if (diag == Diagonal.Unit) dot += x[offsetX + i * incX];
				// else the diagonal is ignored altogether (we operate only on the superdiagonal part)

				// Take into account the rest and store it
				int idxY = offsetY + i * incY;
				y[idxY] = alpha * dot + beta * y[idxY];
			}
		}
	}
}
