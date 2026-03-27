namespace MGroup.LinearAlgebra.Implementations.Managed
{
	using System;
	using System.Diagnostics;

	public partial class ManagedBlasProvider : IBlasProvider
	{
		public void Dgemv(TransposeMatrix trans, int m, int n, double alpha, double[] A, int offsetA, int lda, double[] x, int offsetX, int incX, double beta, double[] y, int offsetY, int incY)
		{
			Debug.Assert(A != null);
			Debug.Assert(x != null);
			Debug.Assert(y != null);
			Debug.Assert(offsetA >= 0);
			Debug.Assert(offsetX >= 0);
			Debug.Assert(offsetY >= 0);
			Debug.Assert(lda >= Math.Max(1, m));
			Debug.Assert(incX != 0);
			Debug.Assert(incY != 0);

			if (m == 0 || n == 0 || (alpha == 0.0 && beta == 1.0))
			{
				return;
			}

			int lenX = (trans == TransposeMatrix.NoTranspose) ? n : m;
			int lenY = (trans == TransposeMatrix.NoTranspose) ? m : n;

			int startX = (incX > 0) ? offsetX : offsetX + (lenX - 1) * incX;
			int startY = (incY > 0) ? offsetY : offsetY + (lenY - 1) * incY;

			if (beta != 1.0)
			{
				int iy = startY;

				if (incY == 1)
				{
					if (beta == 0.0)
					{
						for (int i = 0; i < lenY; i++)
						{
							y[offsetY + i] = 0.0;
						}
					}
					else
					{
						for (int i = 0; i < lenY; i++)
						{
							y[offsetY + i] *= beta;
						}
					}
				}
				else
				{
					if (beta == 0.0)
					{
						for (int i = 0; i < lenY; i++)
						{
							y[iy] = 0.0;
							iy += incY;
						}
					}
					else
					{
						for (int i = 0; i < lenY; i++)
						{
							y[iy] *= beta;
							iy += incY;
						}
					}
				}
			}

			if (alpha == 0.0)
			{
				return;
			}

			if (trans == TransposeMatrix.NoTranspose)
			{
				int ix = startX;

				for (int j = 0; j < n; j++)
				{
					double xj = x[ix];

					if (xj != 0.0)
					{
						double temp = alpha * xj;
						int aColOffset = offsetA + j * lda;
						int iy = startY;

						if (incY == 1)
						{
							for (int i = 0; i < m; i++)
							{
								y[offsetY + i] += temp * A[aColOffset + i];
							}
						}
						else
						{
							for (int i = 0; i < m; i++)
							{
								y[iy] += temp * A[aColOffset + i];
								iy += incY;
							}
						}
					}

					ix += incX;
				}
			}
			else // Transpose
			{
				int iy = startY;

				for (int j = 0; j < n; j++)
				{
					double sum = 0.0;
					int aColOffset = offsetA + j * lda;
					int ix = startX;

					if (incX == 1)
					{
						for (int i = 0; i < m; i++)
						{
							sum += A[aColOffset + i] * x[offsetX + i];
						}
					}
					else
					{
						for (int i = 0; i < m; i++)
						{
							sum += A[aColOffset + i] * x[ix];
							ix += incX;
						}
					}

					y[iy] += alpha * sum;
					iy += incY;
				}
			}
		}

		public void DgemvRowMajor(TransposeMatrix transA, int m, int n, double[] a, double[] x, double[] y)
		{
			if (transA == TransposeMatrix.NoTranspose)
			{
				for (var i = 0; i < m; ++i)
				{
					var rowStart = i * n;
					double sum = 0;
					for (var j = 0; j < n; ++j)
					{
						sum += a[rowStart + j] * x[j];
					}
					y[i] = sum;
				}
			}
			else
			{
				for (var j = 0; j < n; ++j)
				{
					double sum = 0;
					for (var i = 0; i < m; ++i)
					{
						sum += a[i * n + j] * x[i];
					}
					y[j] = sum;
				}
			}
		}

	}
}
