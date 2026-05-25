namespace MGroup.LinearAlgebra.Implementations.Managed
{
	using System;

	using MGroup.LinearAlgebra.Implementations.Managed.Custom;
	//using DotNumerics.LinearAlgebra.CSLapack;
	using MGroup.LinearAlgebra.Implementations.DotNumerics;

	/// <summary>
	/// Provides managed C# implementations of the linear algebra operations defined by <see cref="IBlasProvider"/>. Uses the 
	/// library DotNumerics (see http://www.dotnumerics.com/NumericalLibraries/LinearAlgebra/CSBlas/Default.aspx) for the most
	/// part. For BLAS subroutines not provided by DotNumerics, custom C# implementations are used instead. 
	/// </summary>
	public class DotNumericsBlasProvider : IBlasProvider
	{
		//TODO: perhaps these should not be static.
		private static readonly DAXPY daxpy = new DAXPY();
		private static readonly DCOPY dcopy = new DCOPY();
		private static readonly DDOT ddot = new DDOT();
		private static readonly DGEMM dgemm = new DGEMM();
		private static readonly DGEMV dgemv = new DGEMV();
		private static readonly DNRM2 dnrm2 = new DNRM2();
		private static readonly DSCAL dscal = new DSCAL();
		private static readonly DSWAP dswap = new DSWAP();
		private static readonly DTRSM dtrsm = new DTRSM();
		private static readonly DTRSV dtrsv = new DTRSV();

		private static readonly CustomBlasProvider customBlas = CustomBlasProvider.UniqueInstance;

		public static DotNumericsBlasProvider UniqueInstance { get; } = new DotNumericsBlasProvider();

		private DotNumericsBlasProvider() { } // private constructor for singleton pattern

		#region BLAS Level 1
		public void Daxpby(int n, double alpha, double[] x, int offsetX, int incX, double beta, double[] y, int offsetY, int incY)
		{
			dscal.Run(n, beta, ref y, offsetY, incY);
			daxpy.Run(n, alpha, x, offsetX, incX, ref y, offsetY, incY);
		}

		/// <summary>
		/// See http://www.dotnumerics.com/NumericalLibraries/LinearAlgebra/CSharpCodeFiles/daxpy.aspx
		/// </summary>
		public void Daxpy(int n, double alpha, double[] x, int offsetX, int incX, double[] y, int offsetY, int incY)
			=> daxpy.Run(n, alpha, x, offsetX, incX, ref y, offsetY, incY);

		/// <summary>
		/// See http://www.dotnumerics.com/NumericalLibraries/LinearAlgebra/CSharpCodeFiles/dcopy.aspx
		/// </summary>
		public void Dcopy(int n, double[] x, int offsetX, int incX, double[] y, int offsetY, int incY)
			=> dcopy.Run(n, x, offsetX, incX, ref y, offsetY, incY);

		/// <summary>
		/// See http://www.dotnumerics.com/NumericalLibraries/LinearAlgebra/CSharpCodeFiles/ddot.aspx
		/// </summary>
		public double Ddot(int n, double[] x, int offsetX, int incX, double[] y, int offsetY, int incY)
			=> ddot.Run(n, x, offsetX, incX, y, offsetY, incY);

		/// <summary>
		/// See http://www.dotnumerics.com/NumericalLibraries/LinearAlgebra/CSharpCodeFiles/dnrm2.aspx
		/// </summary>
		public double Dnrm2(int n, double[] x, int offsetX, int incX)
			=> dnrm2.Run(n, x, offsetX, incX);

		/// <summary>
		/// See http://www.dotnumerics.com/NumericalLibraries/LinearAlgebra/CSharpCodeFiles/dscal.aspx
		/// </summary>
		public void Dscal(int n, double alpha, double[] x, int offsetX, int incX)
			=> dscal.Run(n, alpha, ref x, offsetX, incX);

		/// <summary>
		/// See http://www.dotnumerics.com/NumericalLibraries/LinearAlgebra/CSharpCodeFiles/dswap.aspx
		/// </summary>
		public void Dswap(int n, double[] x, int offsetX, int incX, double[] y, int offsetY, int incY) 
			=> dswap.Run(n, ref x, offsetX, incX, ref y, offsetY, incY);

		#endregion

		#region BLAS Level 2

		/// <summary>
		/// See http://www.dotnumerics.com/NumericalLibraries/LinearAlgebra/CSharpCodeFiles/dgemv.aspx
		/// </summary>
		public void Dgemv(TransposeMatrix transA, int m, int n,
			double alpha, double[] a, int offsetA, int ldA, double[] x, int offsetX, int incX,
			double beta, double[] y, int offsetY, int incY)
			=> dgemv.Run(transA.Translate(), m, n, alpha, a, offsetA, ldA, x, offsetX, incX, beta, ref y, offsetY, incY);

		public void DgemvRowMajor(TransposeMatrix transA, int m, int n, double[] a, double[] x, double[] y)
		{
			customBlas.DgemvRowMajor(transA, m, n, a, x, y);
		}

		public void Dspmv(StoredTriangle uplo, int n,
			double alpha, double[] a, int offsetA, double[] x, int offsetX, int incX,
			double beta, double[] y, int offsetY, int incY)
		{
			customBlas.Dspmv(uplo, n, alpha, a, offsetA, x, offsetX, incX, beta, y, offsetY, incY);
		}

		public void Dtpmv(StoredTriangle uplo, TransposeMatrix transA, DiagonalValues diag, int n,
			double[] a, int offsetA, double[] x, int offsetX, int incX)
		{
			customBlas.Dtpmv(uplo, transA, diag, n, a, offsetA, x, offsetX, incX);
		}

		public void Dtpsv(StoredTriangle uplo, TransposeMatrix transA, DiagonalValues diag, int n,
			double[] a, int offsetA, double[] x, int offsetX, int incX)
		{
			customBlas.Dtpsv(uplo, transA, diag, n, a, offsetA, x, offsetX, incX);
		}

		/// <summary>
		/// See http://www.dotnumerics.com/NumericalLibraries/LinearAlgebra/CSharpCodeFiles/dtrsv.aspx
		/// </summary>
		public void Dtrsv(StoredTriangle uplo, TransposeMatrix transA, DiagonalValues diag, int n,
			double[] a, int offsetA, int ldA, double[] x, int offsetX, int incX)
			=> dtrsv.Run(uplo.Translate(), transA.Translate(), diag.Translate(), n, a, offsetA, ldA, ref x, offsetX, incX);
		#endregion

		#region BLAS Level 3

		/// <summary>
		/// See http://www.dotnumerics.com/NumericalLibraries/LinearAlgebra/CSharpCodeFiles/dgemm.aspx
		/// </summary>
		public void Dgemm(TransposeMatrix transA, TransposeMatrix transB, int m, int n, int k, double alpha,
			double[] a, int offsetA, int ldA, double[] b, int offsetB, int ldB, double beta, double[] c, int offsetC, int ldC)
			=> dgemm.Run(transA.Translate(), transB.Translate(), m, n, k, alpha, a, offsetA, ldA, b, offsetB, ldB,
				beta, ref c, offsetC, ldC);

		/// <summary>
		/// See http://www.dotnumerics.com/NumericalLibraries/LinearAlgebra/CSharpCodeFiles/dtrsm.aspx
		/// </summary>
		public void Dtrsm(MultiplicationSide side, StoredTriangle uplo, TransposeMatrix transA, DiagonalValues diag, int m, int n, double alpha, double[] a, int offsetA, int ldA, double[] b, int offsetB, int ldB) 
			=> dtrsm.Run(side.Translate(), uplo.Translate(), transA.Translate(), diag.Translate(),
				m, n, alpha, a, offsetA, ldA, ref b, offsetB, ldB);
		#endregion
	}
}
