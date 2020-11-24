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
		/// The imaginary parts of eigenvalues of a matrix. Will be empty if all eigenvalues are real, 
		/// e.g. for symmetric matrices.
		/// </summary>
		Vector EigenvaluesImaginary { get; }

		/// <summary>
		/// The right eigenvectors of a matrix. Each column j is an eigenvector v that corresponds to 
		/// lambda=<see cref="Eigenvalues"/>[j], such as A*v = lambda * v.
		/// This matrix may be empty, if the user requested that no right eigenvectors were computed.
		/// </summary>
		Matrix EigenvectorsRight { get; }

		/// <summary>
		/// The left eigenvectors of a matrix. Each column j is an eigenvector v that corresponds to 
		/// lambda=<see cref="Eigenvalues"/>[j], such as transpose(v)*A = lambda * v.
		/// This matrix may be empty, if the user requested that no left eigenvectors were computed.
		/// </summary>
		Matrix EigenvectorsLeft { get; }

		/// <summary>
		/// The number of rows and columns of the original (square) matrix.
		/// </summary>
		int Order { get; }
	}
}
