namespace MGroup.LinearAlgebra.Implementations.Managed.Custom
{
	using System;
	using System.Diagnostics;

	public partial class CustomBlasProvider : IBlasProvider
	{
		public void Dtrsm(MultiplicationSide side, StoredTriangle uplo, TransposeMatrix transA, DiagonalValues diag, int m, int n, double alpha, double[] a, int offsetA, int ldA, double[] b, int offsetB, int ldB)
		{
			Debug.Assert(a != null);
			Debug.Assert(b != null);
			Debug.Assert(offsetA >= 0);
			Debug.Assert(offsetB >= 0);

			bool leftSide = (side == MultiplicationSide.Left);
			bool upper = (uplo == StoredTriangle.Upper);
			bool noTrans = (transA == TransposeMatrix.NoTranspose);
			bool nonUnit = (diag == DiagonalValues.NonUnit);

			int nrowA = leftSide ? m : n;
			Debug.Assert(ldA >= Math.Max(1, nrowA));
			Debug.Assert(ldB >= Math.Max(1, m));

			if (m == 0 || n == 0)
			{
				return;
			}

			if (alpha == 0.0)
			{
				for (int j = 0; j < n; j++)
				{
					int colB = offsetB + j * ldB;
					for (int i = 0; i < m; i++)
					{
						b[colB + i] = 0.0;
					}
				}
				return;
			}

			if (leftSide)
			{
				if (noTrans)
				{
					if (upper)
					{
						for (int j = 0; j < n; j++)
						{
							int colB = offsetB + j * ldB;

							if (alpha != 1.0)
							{
								for (int i = 0; i < m; i++) b[colB + i] *= alpha;
							}

							for (int k = m - 1; k >= 0; k--)
							{
								int idx = colB + k;
								if (b[idx] != 0.0)
								{
									if (nonUnit)
									{
										b[idx] /= a[offsetA + k * ldA + k];
									}

									double temp = b[idx];
									int colA = offsetA + k * ldA;

									for (int i = 0; i < k; i++)
									{
										b[colB + i] -= temp * a[colA + i];
									}
								}
							}
						}
					}
					else
					{
						for (int j = 0; j < n; j++)
						{
							int colB = offsetB + j * ldB;

							if (alpha != 1.0)
							{
								for (int i = 0; i < m; i++) b[colB + i] *= alpha;
							}

							for (int k = 0; k < m; k++)
							{
								int idx = colB + k;
								if (b[idx] != 0.0)
								{
									if (nonUnit)
									{
										b[idx] /= a[offsetA + k * ldA + k];
									}

									double temp = b[idx];
									int colA = offsetA + k * ldA;

									for (int i = k + 1; i < m; i++)
									{
										b[colB + i] -= temp * a[colA + i];
									}
								}
							}
						}
					}
				}
				else
				{
					if (upper)
					{
						for (int j = 0; j < n; j++)
						{
							int colB = offsetB + j * ldB;

							for (int i = 0; i < m; i++)
							{
								double sum = alpha * b[colB + i];
								int rowA = offsetA + i * ldA;

								for (int k = 0; k < i; k++)
								{
									sum -= a[rowA + k] * b[colB + k];
								}

								if (nonUnit)
								{
									sum /= a[rowA + i];
								}

								b[colB + i] = sum;
							}
						}
					}
					else
					{
						for (int j = 0; j < n; j++)
						{
							int colB = offsetB + j * ldB;

							for (int i = m - 1; i >= 0; i--)
							{
								double sum = alpha * b[colB + i];
								int rowA = offsetA + i * ldA;

								for (int k = i + 1; k < m; k++)
								{
									sum -= a[rowA + k] * b[colB + k];
								}

								if (nonUnit)
								{
									sum /= a[rowA + i];
								}

								b[colB + i] = sum;
							}
						}
					}
				}
			}
			else
			{
				if (noTrans)
				{
					if (upper)
					{
						for (int j = 0; j < n; j++)
						{
							int colB = offsetB + j * ldB;

							if (alpha != 1.0)
							{
								for (int i = 0; i < m; i++) b[colB + i] *= alpha;
							}

							for (int k = 0; k < j; k++)
							{
								double aVal = a[offsetA + j * ldA + k];
								if (aVal != 0.0)
								{
									int colBk = offsetB + k * ldB;
									for (int i = 0; i < m; i++)
									{
										b[colB + i] -= aVal * b[colBk + i];
									}
								}
							}

							if (nonUnit)
							{
								double inv = 1.0 / a[offsetA + j * ldA + j];
								for (int i = 0; i < m; i++)
								{
									b[colB + i] *= inv;
								}
							}
						}
					}
					else
					{
						for (int j = n - 1; j >= 0; j--)
						{
							int colB = offsetB + j * ldB;

							if (alpha != 1.0)
							{
								for (int i = 0; i < m; i++) b[colB + i] *= alpha;
							}

							for (int k = j + 1; k < n; k++)
							{
								double aVal = a[offsetA + j * ldA + k];
								if (aVal != 0.0)
								{
									int colBk = offsetB + k * ldB;
									for (int i = 0; i < m; i++)
									{
										b[colB + i] -= aVal * b[colBk + i];
									}
								}
							}

							if (nonUnit)
							{
								double inv = 1.0 / a[offsetA + j * ldA + j];
								for (int i = 0; i < m; i++)
								{
									b[colB + i] *= inv;
								}
							}
						}
					}
				}
				else
				{
					if (upper)
					{
						for (int k = n - 1; k >= 0; k--)
						{
							if (nonUnit)
							{
								double inv = 1.0 / a[offsetA + k * ldA + k];
								int colBk = offsetB + k * ldB;
								for (int i = 0; i < m; i++)
								{
									b[colBk + i] *= inv;
								}
							}

							for (int j = 0; j < k; j++)
							{
								double aVal = a[offsetA + k * ldA + j];
								if (aVal != 0.0)
								{
									int colBj = offsetB + j * ldB;
									int colBk = offsetB + k * ldB;
									for (int i = 0; i < m; i++)
									{
										b[colBj + i] -= aVal * b[colBk + i];
									}
								}
							}

							if (alpha != 1.0)
							{
								int colBk = offsetB + k * ldB;
								for (int i = 0; i < m; i++)
								{
									b[colBk + i] *= alpha;
								}
							}
						}
					}
					else
					{
						for (int k = 0; k < n; k++)
						{
							if (nonUnit)
							{
								double inv = 1.0 / a[offsetA + k * ldA + k];
								int colBk = offsetB + k * ldB;
								for (int i = 0; i < m; i++)
								{
									b[colBk + i] *= inv;
								}
							}

							for (int j = k + 1; j < n; j++)
							{
								double aVal = a[offsetA + k * ldA + j];
								if (aVal != 0.0)
								{
									int colBj = offsetB + j * ldB;
									int colBk = offsetB + k * ldB;
									for (int i = 0; i < m; i++)
									{
										b[colBj + i] -= aVal * b[colBk + i];
									}
								}
							}

							if (alpha != 1.0)
							{
								int colBk = offsetB + k * ldB;
								for (int i = 0; i < m; i++)
								{
									b[colBk + i] *= alpha;
								}
							}
						}
					}
				}
			}
		}
	}
}
