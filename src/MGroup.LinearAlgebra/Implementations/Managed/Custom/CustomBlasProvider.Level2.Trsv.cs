namespace MGroup.LinearAlgebra.Implementations.Managed.Custom
{
	using System;
	using System.Diagnostics;

	public partial class CustomBlasProvider : IBlasProvider
	{
		public void Dtrsv(StoredTriangle uplo, TransposeMatrix trans, DiagonalValues diag, int n, double[] A, int offsetA, int lda, double[] x, int offsetX, int incX)
		{
			Debug.Assert(A != null);
			Debug.Assert(x != null);
			Debug.Assert(offsetA >= 0);
			Debug.Assert(offsetX >= 0);
			Debug.Assert(lda >= Math.Max(1, n));
			Debug.Assert(incX != 0);

			if (n == 0)
			{
				return;
			}

			bool nonUnit = (diag == DiagonalValues.NonUnit);

			int startX;
			if (incX <= 0)
			{
				startX = offsetX + (n - 1) * incX;
			}
			else
			{
				startX = offsetX;
			}

			if (trans == TransposeMatrix.NoTranspose)
			{
				if (uplo == StoredTriangle.Upper)
				{
					if (incX == 1)
					{
						for (int j = n - 1; j >= 0; j--)
						{
							int xIndex = offsetX + j;

							if (x[xIndex] != 0.0)
							{
								if (nonUnit)
								{
									x[xIndex] /= A[offsetA + j * lda + j];
								}

								double temp = x[xIndex];
								int colOffset = offsetA + j * lda;

								for (int i = j - 1; i >= 0; i--)
								{
									x[offsetX + i] -= temp * A[colOffset + i];
								}
							}
						}
					}
					else
					{
						int jx = startX + (n - 1) * incX;

						for (int j = n - 1; j >= 0; j--)
						{
							if (x[jx] != 0.0)
							{
								if (nonUnit)
								{
									x[jx] /= A[offsetA + j * lda + j];
								}

								double temp = x[jx];
								int ix = jx;
								int colOffset = offsetA + j * lda;

								for (int i = j - 1; i >= 0; i--)
								{
									ix -= incX;
									x[ix] -= temp * A[colOffset + i];
								}
							}

							jx -= incX;
						}
					}
				}
				else // Lower
				{
					if (incX == 1)
					{
						for (int j = 0; j < n; j++)
						{
							int xIndex = offsetX + j;

							if (x[xIndex] != 0.0)
							{
								if (nonUnit)
								{
									x[xIndex] /= A[offsetA + j * lda + j];
								}

								double temp = x[xIndex];
								int colOffset = offsetA + j * lda;

								for (int i = j + 1; i < n; i++)
								{
									x[offsetX + i] -= temp * A[colOffset + i];
								}
							}
						}
					}
					else
					{
						int jx = startX;

						for (int j = 0; j < n; j++)
						{
							if (x[jx] != 0.0)
							{
								if (nonUnit)
								{
									x[jx] /= A[offsetA + j * lda + j];
								}

								double temp = x[jx];
								int ix = jx;
								int colOffset = offsetA + j * lda;

								for (int i = j + 1; i < n; i++)
								{
									ix += incX;
									x[ix] -= temp * A[colOffset + i];
								}
							}

							jx += incX;
						}
					}
				}
			}
			else // Transpose
			{
				if (uplo == StoredTriangle.Upper)
				{
					if (incX == 1)
					{
						for (int j = 0; j < n; j++)
						{
							double temp = x[offsetX + j];
							int colOffset = offsetA + j * lda;

							for (int i = 0; i < j; i++)
							{
								temp -= A[colOffset + i] * x[offsetX + i];
							}

							if (nonUnit)
							{
								temp /= A[offsetA + j * lda + j];
							}

							x[offsetX + j] = temp;
						}
					}
					else
					{
						int jx = startX;

						for (int j = 0; j < n; j++)
						{
							double temp = x[jx];
							int ix = startX;
							int colOffset = offsetA + j * lda;

							for (int i = 0; i < j; i++)
							{
								temp -= A[colOffset + i] * x[ix];
								ix += incX;
							}

							if (nonUnit)
							{
								temp /= A[offsetA + j * lda + j];
							}

							x[jx] = temp;
							jx += incX;
						}
					}
				}
				else // Lower
				{
					if (incX == 1)
					{
						for (int j = n - 1; j >= 0; j--)
						{
							double temp = x[offsetX + j];
							int colOffset = offsetA + j * lda;

							for (int i = n - 1; i > j; i--)
							{
								temp -= A[colOffset + i] * x[offsetX + i];
							}

							if (nonUnit)
							{
								temp /= A[offsetA + j * lda + j];
							}

							x[offsetX + j] = temp;
						}
					}
					else
					{
						int start = startX + (n - 1) * incX;
						int jx = start;

						for (int j = n - 1; j >= 0; j--)
						{
							double temp = x[jx];
							int ix = start;
							int colOffset = offsetA + j * lda;

							for (int i = n - 1; i > j; i--)
							{
								temp -= A[colOffset + i] * x[ix];
								ix -= incX;
							}

							if (nonUnit)
							{
								temp /= A[offsetA + j * lda + j];
							}

							x[jx] = temp;
							jx -= incX;
						}
					}
				}
			}
		}
	}
}
