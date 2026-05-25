namespace MGroup.LinearAlgebra.Implementations.Managed.Custom
{
	using System;

	using MGroup.LinearAlgebra.Commons;

	/// <summary>
	/// Provides custom managed C# implementations of the linear algebra operations defined by <see cref="ILapackProvider"/>.  
	/// </summary>
	public partial class CustomLapackProvider : ILapackProvider
	{
		private CustomBlasProvider blas = CustomBlasProvider.UniqueInstance;

		public static CustomLapackProvider UniqueInstance { get; } = new CustomLapackProvider();

		private CustomLapackProvider() { } // private constructor for singleton pattern

		public void Dpotrf(string uplo, int n, double[] a, int offsetA, int ldA, ref int info)
		{
			try
			{
				if (IsUpper(uplo)) LapackImplementations.CholeskyUpperFullColMajor(n, a, offsetA, ldA, ref info);
				else LapackImplementations.CholeskyLowerFullColMajor(n, a, offsetA, ldA, ref info);
			}
			catch (ArgumentException)
			{
				info = -1;
			}
		}

		public void Dpotri(string uplo, int n, double[] a, int offsetA, int ldA, ref int info)
		{
			if (ldA != n) throw new NotImplementedException("dpotri only works for ldA=n");

			// Start with an identity matrix
			var inverse = new double[n * n];
			for (int i = 0; i < n; ++i) inverse[i * ldA + i] = 1.0;

			// Solve (L*L^T) * inverse = I or (U^T*U) * inverse = I
			int infoSolve = LapackUtilities.DefaultInfo;
			Dpotrs(uplo, n, n, a, offsetA, ldA, inverse, 0, n, ref infoSolve);

			// Copy the inverse matrix over the factorization
			Array.Copy(inverse, 0, a, offsetA, n * n);

			info = 0; //TODO add checks
		}

		public void Dpotrs(string uplo, int n, int nRhs, double[] a, int offsetA, int ldA,
			double[] b, int offsetB, int ldB, ref int info)
		{
			try
			{
				if (IsUpper(uplo)) // A*X=B <=> U^T * (U*X) = B
				{
					blas.Dtrsm(MultiplicationSide.Left, StoredTriangle.Upper, TransposeMatrix.Transpose, DiagonalValues.NonUnit, n, nRhs, 1.0, a, offsetA, ldA, b, offsetB, ldB); // B = U^T \ B
					blas.Dtrsm(MultiplicationSide.Left, StoredTriangle.Upper, TransposeMatrix.NoTranspose, DiagonalValues.NonUnit, n, nRhs, 1.0, a, offsetA, ldA, b, offsetB, ldB); // B = U \ B
				}
				else // A*X=B <=> L * (L^T*X) = B
				{
					blas.Dtrsm(MultiplicationSide.Left, StoredTriangle.Lower, TransposeMatrix.NoTranspose, DiagonalValues.NonUnit, n, nRhs, 1.0, a, offsetA, ldA, b, offsetB, ldB); // B = L \ B
					blas.Dtrsm(MultiplicationSide.Left, StoredTriangle.Lower, TransposeMatrix.Transpose, DiagonalValues.NonUnit, n, nRhs, 1.0, a, offsetA, ldA, b, offsetB, ldB); // B = L^T \ B
				}
				info = 0; // TODO: needs more checks
			}
			catch (ArgumentException)
			{
				info = -1;
			}
		}

		public void Dpptrf(string uplo, int n, double[] a, int offsetA, ref int info)
		{
			try
			{
				if (IsUpper(uplo)) LapackImplementations.CholeskyUpperPackedColMajor(n, a, offsetA, ref info);
				else LapackImplementations.CholeskyLowerPackedColMajor(n, a, offsetA, ref info);
			}
			catch (ArgumentException)
			{
				info = -1;
			}
		}

		public void Dpptri(string uplo, int n, double[] a, int offsetA, ref int info)
		{
			try
			{
				// Start with an identity matrix
				var inverse = new double[n * n];
				for (int i = 0; i < n; ++i) inverse[i * n + i] = 1.0;

				// Solve (L*L^T) * inverse = I or (U^T*U) * inverse = I
				int infoSolve = LapackUtilities.DefaultInfo;
				Dpptrs(uplo, n, n, a, offsetA, inverse, 0, n, ref infoSolve);

				// Copy the inverse matrix over the factorization
				if (IsUpper(uplo)) Conversions.FullColMajorToPackedUpperColMajor(n, inverse, a, offsetA);
				else Conversions.FullColMajorToPackedLowerColMajor(n, inverse, a, offsetA);
				info = 0; // TODO: needs more checks
			}
			catch (ArgumentException)
			{
				info = -1;
			}
		}

		public void Dpptrs(string uplo, int n, int nRhs, double[] a, int offsetA, double[] b, int offsetB, int ldB, ref int info)
		{
			try
			{
				if (IsUpper(uplo))
				{
					// Process each column separately
					for (int i = 0; i < nRhs; ++i)
					{
						// b = U^T \ b
						CustomBlasProvider.UniqueInstance.Dtpsv(StoredTriangle.Upper, TransposeMatrix.Transpose,
							DiagonalValues.NonUnit, n, a, offsetA, b, offsetB + i * nRhs, 1);

						// b = U \ b
						CustomBlasProvider.UniqueInstance.Dtpsv(StoredTriangle.Upper, TransposeMatrix.NoTranspose,
							DiagonalValues.NonUnit, n, a, offsetA, b, offsetB + i * nRhs, 1);
					}
				}
				else
				{
					// Process each column separately
					for (int i = 0; i < nRhs; ++i)
					{
						// b = L \ b
						CustomBlasProvider.UniqueInstance.Dtpsv(StoredTriangle.Lower, TransposeMatrix.NoTranspose,
							DiagonalValues.NonUnit, n, a, offsetA, b, offsetB + i * n, 1);

						// b = L^T \ b
						CustomBlasProvider.UniqueInstance.Dtpsv(StoredTriangle.Lower, TransposeMatrix.Transpose,
							DiagonalValues.NonUnit, n, a, offsetA, b, offsetB + i * n, 1);
					}
				}
				info = 0; // TODO: needs more checks
			}
			catch (ArgumentException)
			{
				info = -1;
			}
		}

		private static bool IsUpper(string uplo)
		{
			if (uplo.Equals("L") || uplo.Equals("l")) return false;
			else if (uplo.Equals("U") || uplo.Equals("u")) return true;
			else throw new ArgumentException("Parameter uplo must be U, u, L or l");
		}

		public void Dgeev(string jobVl, string jobVr, int n, ref double[] a, int offsetA, int ldA, ref double[] wr, int offsetWr, ref double[] wi, int offsetWi, ref double[] vl, int offsetVl, int ldVl, ref double[] vr, int offsetVr, int ldVr, ref double[] work, int offsetWork, int lWork, ref int info) => throw new NotImplementedException();
		
		public void Dgelqf(int m, int n, double[] a, int offsetA, int ldA, double[] tau, int offsetTau, double[] work, int offsetWork, int lWork, ref int info) => throw new NotImplementedException();
		
		public void Dgeqrf(int m, int n, double[] a, int offsetA, int ldA, double[] tau, int offsetTau, double[] work, int offsetWork, int lWork, ref int info) => throw new NotImplementedException();
		
		public void Dgetri(int n, double[] a, int offsetA, int ldA, int[] ipiv, int offsetIpiv, double[] work, int offsetWork, int lWork, ref int info) => throw new NotImplementedException();
		
		public void Dgetrs(string transA, int n, int nRhs, double[] a, int offsetA, int ldA, int[] ipiv, int offsetIpiv, double[] b, int offsetB, int ldB, ref int info) => throw new NotImplementedException();
		
		public void Dorglq(int m, int n, int k, double[] a, int offsetA, int ldA, double[] tau, int offsetTau, double[] work, int offsetWork, int lWork, ref int info) => throw new NotImplementedException();
		
		public void Dorgqr(int m, int n, int k, double[] a, int offsetA, int ldA, double[] tau, int offsetTau, double[] work, int offsetWork, int lWork, ref int info) => throw new NotImplementedException();
		
		public void Dormlq(string side, string transQ, int m, int n, int k, double[] a, int offsetA, int ldA, double[] tau, int offsetTau, double[] c, int offsetC, int ldC, double[] work, int offsetWork, int lWork, ref int info) => throw new NotImplementedException();
		
		public void Dormqr(string side, string transQ, int m, int n, int k, double[] a, int offsetA, int ldA, double[] tau, int offsetTau, double[] c, int offsetC, int ldC, double[] work, int offsetWork, int lWork, ref int info) => throw new NotImplementedException();

		public void Dsyev(string jobz, string uplo, int n, ref double[] a, int offsetA, int ldA, ref double[] w, int offsetW, ref double[] work, int offsetWork, int lWork, ref int info) => throw new NotImplementedException();
	}
}
