namespace MGroup.LinearAlgebra.Implementations.Managed.Custom
{
	using System;
	using System.Diagnostics;

	public partial class CustomBlasProvider : IBlasProvider
	{
		public void Dgemm(TransposeMatrix transA, TransposeMatrix transB, int m, int n, int k, double alpha, double[] A, int offsetA, int lda, double[] B, int offsetB, int ldb, double beta, double[] C, int offsetC, int ldc)
		{
			Debug.Assert(A != null);
			Debug.Assert(B != null);
			Debug.Assert(C != null);
			Debug.Assert(offsetA >= 0);
			Debug.Assert(offsetB >= 0);
			Debug.Assert(offsetC >= 0);

			bool notransA = (transA == TransposeMatrix.NoTranspose);
			bool notransB = (transB == TransposeMatrix.NoTranspose);

			int nrowA = notransA ? m : k;
			int nrowB = notransB ? k : n;

			Debug.Assert(lda >= Math.Max(1, nrowA));
			Debug.Assert(ldb >= Math.Max(1, nrowB));
			Debug.Assert(ldc >= Math.Max(1, m));

			if (m == 0 || n == 0 || ((alpha == 0.0 || k == 0) && beta == 1.0))
			{
				return;
			}

			if (alpha == 0.0)
			{
				if (beta == 0.0)
				{
					for (int j = 0; j < n; j++)
					{
						int cCol = offsetC + j * ldc;
						for (int i = 0; i < m; i++)
						{
							C[cCol + i] = 0.0;
						}
					}
				}
				else
				{
					for (int j = 0; j < n; j++)
					{
						int cCol = offsetC + j * ldc;
						for (int i = 0; i < m; i++)
						{
							C[cCol + i] *= beta;
						}
					}
				}
				return;
			}

			if (notransB)
			{
				if (notransA)
				{
					for (int j = 0; j < n; j++)
					{
						int cCol = offsetC + j * ldc;

						if (beta == 0.0)
						{
							for (int i = 0; i < m; i++) C[cCol + i] = 0.0;
						}
						else if (beta != 1.0)
						{
							for (int i = 0; i < m; i++) C[cCol + i] *= beta;
						}

						for (int l = 0; l < k; l++)
						{
							double bVal = B[offsetB + j * ldb + l];
							if (bVal != 0.0)
							{
								double temp = alpha * bVal;
								int aCol = offsetA + l * lda;

								for (int i = 0; i < m; i++)
								{
									C[cCol + i] += temp * A[aCol + i];
								}
							}
						}
					}
				}
				else
				{
					for (int j = 0; j < n; j++)
					{
						for (int i = 0; i < m; i++)
						{
							double sum = 0.0;
							int aRow = offsetA + i * lda;
							int bCol = offsetB + j * ldb;

							for (int l = 0; l < k; l++)
							{
								sum += A[aRow + l] * B[bCol + l];
							}

							int cIndex = offsetC + j * ldc + i;

							if (beta == 0.0)
							{
								C[cIndex] = alpha * sum;
							}
							else
							{
								C[cIndex] = alpha * sum + beta * C[cIndex];
							}
						}
					}
				}
			}
			else
			{
				if (notransA)
				{
					for (int j = 0; j < n; j++)
					{
						int cCol = offsetC + j * ldc;

						if (beta == 0.0)
						{
							for (int i = 0; i < m; i++) C[cCol + i] = 0.0;
						}
						else if (beta != 1.0)
						{
							for (int i = 0; i < m; i++) C[cCol + i] *= beta;
						}

						for (int l = 0; l < k; l++)
						{
							double bVal = B[offsetB + l * ldb + j];
							if (bVal != 0.0)
							{
								double temp = alpha * bVal;
								int aCol = offsetA + l * lda;

								for (int i = 0; i < m; i++)
								{
									C[cCol + i] += temp * A[aCol + i];
								}
							}
						}
					}
				}
				else
				{
					for (int j = 0; j < n; j++)
					{
						for (int i = 0; i < m; i++)
						{
							double sum = 0.0;
							int aRow = offsetA + i * lda;

							for (int l = 0; l < k; l++)
							{
								sum += A[aRow + l] * B[offsetB + l * ldb + j];
							}

							int cIndex = offsetC + j * ldc + i;

							if (beta == 0.0)
							{
								C[cIndex] = alpha * sum;
							}
							else
							{
								C[cIndex] = alpha * sum + beta * C[cIndex];
							}
						}
					}
				}
			}
		}
	}
}
