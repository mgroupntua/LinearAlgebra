using System;
using System.Collections.Generic;
using System.Text;
using MGroup.LinearAlgebra.Exceptions;
using static MGroup.LinearAlgebra.Providers.LapackUtilities;

namespace MGroup.LinearAlgebra.Providers
{
	/// <summary>
	/// Simplifies the use of LAPACK (see <see cref="ILapackProvider"/>) linear algebra operations that concern the calculation
	/// of eigenvalues, eigenvectors, singular values, etc. with double precision arithmetic. Such simplifications are error 
	/// checking, handling workspace arrays, enums instead of string arguments etc. This facade is meant to provide a managed 
	/// object-oriented alternative the LAPACKE library used in C.
	/// Authors: Serafeim Bakalakos
	/// </summary>
	internal class LapackEigensystemsFacade
	{
		internal LapackEigensystemsFacade(ILapackProvider provider)
		{
			this.Provider = provider;
		}

		internal ILapackProvider Provider { get; }

		internal void Dgeev(EigensystemJob jobLeft, EigensystemJob jobRight, int orderA, 
			double[] matrixA, int offsetA, int leadingDimA, 
			double[] eigvaluesReal, int offsetEigvaluesReal, double[] eigvaluesImaginary, int offsetEigvaluesImaginary,
			double[] eigvectorsLeft, int offsetEigvectorsLeft, int leadingDimEigvectorsLeft,
			double[] eigvectorsRight, int offsetEigvectorsRight, int leadingDimEigvectorsRight)
		{
			int info = DefaultInfo;
			QueryWorkspaceAndExecute((work, offsetWork, lWork) => Provider.Dgeev(jobLeft.Translate(), jobRight.Translate(), 
				orderA, ref matrixA, offsetA, leadingDimA,
				ref eigvaluesReal, offsetEigvaluesReal, ref eigvaluesImaginary, offsetEigvaluesImaginary,
				ref eigvectorsLeft, offsetEigvectorsLeft, leadingDimEigvectorsLeft,
				ref eigvectorsRight, offsetEigvectorsRight, leadingDimEigvectorsRight,
				ref work, offsetWork, lWork, ref info));

			if (info > 0)
			{
				throw new LapackException($"The QR algorithm failed to compute all the eigenvalues, and no eigenvectors have "
				+$"been computed. Elements {info+1:N} of the eigenvalues arrays (real and imaginary) contain eigenvalues which" 
				+ "have converged");
			}
			else if (info < 0) ProcessNegativeInfo(info);
		}

		internal void Dsyev(EigensystemJob job, StoredTriangle triangle, int orderA, double[] matrixA, int offsetA,
			int leadingDimA, double[] eigenvalues, int offsetEigenvalues)
		{
			int info = DefaultInfo;
			QueryWorkspaceAndExecute((work, offsetWork, lWork) => Provider.Dsyev(
				job.Translate(), triangle.Translate(), orderA, ref matrixA, offsetA, leadingDimA, ref eigenvalues,
				offsetEigenvalues, ref work, offsetWork, lWork, ref info));

			if (info > 0)
			{
				throw new LapackException($"The algorithm failed to converge. There were {info} elements of an intermediate"
					+ " tridiagonal form which did not converge to zero");
			}
			else if (info < 0) ProcessNegativeInfo(info);
		}
	}
}
