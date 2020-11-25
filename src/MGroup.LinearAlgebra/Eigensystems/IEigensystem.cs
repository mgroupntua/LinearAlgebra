using System;
using System.Collections.Generic;
using System.Text;
using MGroup.LinearAlgebra.Matrices;
using MGroup.LinearAlgebra.Vectors;

namespace MGroup.LinearAlgebra.Eigensystems
{
	/// <summary>
	/// Represents the eigvalues and eigenvectors of a matrix.
	/// </summary>
	public interface IEigensystem
	{
		/// <summary>
		/// The real parts of eigenvalues of a matrix.
		/// </summary>
		Vector EigenvaluesReal { get; }

		/// <summary>
		/// The imaginary parts of eigenvalues of a matrix. It may be null or equal to zero if all eigenvalues are real, 
		/// e.g. for symmetric matrices.
		/// </summary>
		Vector EigenvaluesImaginary { get; }

		/// <summary>
		/// The left eigenvectors of a matrix. Each column j is an eigenvector v that corresponds to 
		/// lambda=<see cref="Eigenvalues"/>[j], such as transpose(v)*A = lambda * v.
		/// If eigenvalue j is real, then the corresponding left eigenvector is equal to <see cref="EigenvectorsLeft"/>[j].
		/// If eigenvalues j, j+1 are a complex conjugate pair, then the corresponding left eigenvectors are: 
		/// u(j) = <see cref="EigenvectorsLeft"/>[j] + i * <see cref="EigenvectorsLeft"/>[j+1] and
		/// u(j+1) = <see cref="EigenvectorsLeft"/>[j] - i * <see cref="EigenvectorsLeft"/>[j+1].
		/// It may be null, if the user requested that no left eigenvectors were computed.
		/// </summary>
		Matrix EigenvectorsLeft { get; }

		/// <summary>
		/// The right eigenvectors of a matrix. Each column j is an eigenvector v that corresponds to 
		/// lambda=<see cref="Eigenvalues"/>[j], such as A*v = lambda * v.
		/// If eigenvalue j is real, then the corresponding right eigenvector is equal to <see cref="EigenvectorsRight"/>[j].
		/// If eigenvalues j, j+1 are a complex conjugate pair, then the corresponding right eigenvectors are: 
		/// u(j) = <see cref="EigenvectorsRight"/>[j] + i * <see cref="EigenvectorsRight"/>[j+1] and
		/// u(j+1) = <see cref="EigenvectorsRight"/>[j] - i * <see cref="EigenvectorsRight"/>[j+1].
		/// It may be null, if the user requested that no right eigenvectors were computed.
		/// </summary>
		Matrix EigenvectorsRight { get; }

		/// <summary>
		/// The number of rows and columns of the original square matrix.
		/// </summary>
		int Order { get; }
	}
}
