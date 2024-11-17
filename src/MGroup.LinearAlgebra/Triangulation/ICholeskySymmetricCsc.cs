namespace MGroup.LinearAlgebra.Triangulation
{
	using System;
	using System.Collections.Generic;
	using System.Text;
	using MGroup.LinearAlgebra.Exceptions;
	using MGroup.LinearAlgebra.Matrices;

	/// <summary>
	/// Cholesky factorization of a sparse symmetric positive definite matrix. The original matrix must be in Compressed Sparse 
	/// Columns (CSC) format, with only the upper triangle stored.
	/// </summary>
	public interface ICholeskySymmetricCsc : ITriangulation, IDisposable
	{
		/// <summary>
		/// Performs the Cholesky factorization: A = L * L^T of a symmetric positive definite matrix A.
		/// Only the upper triangle of the original matrix is required and is provided in symmetric CSC format by
		/// <paramref name="cscValues"/>, <paramref name="cscRowIndices"/> and <paramref name="cscColOffsets"/>.
		/// </summary>
		/// <param name="order">The number of rows/columns of the square matrix.</param>
		/// <param name="numNonZerosUpper">The number of explicitly stored entries in the upper triangle of the matrix.</param>
		/// <param name="cscValues">
		/// Contains the non-zero entries of the upper triangle. Its length must be equal to <paramref name="numNonZerosUpper"/>.
		/// The non-zero entries of each row must appear consecutively in <paramref name="cscValues"/>. They should also be
		/// sorted in increasing order of their row indices, to speed up subsequent the factorization.
		/// </param>
		/// <param name="cscRowIndices">
		/// Contains the row indices of the non-zero entries. Its length must be equal to <paramref name="numNonZerosUpper"/>.
		/// There is an 1 to 1 matching between these two arrays: <paramref name="cscRowIndices"/>[i] is the row index of the
		/// entry <paramref name="cscValues"/>[i]. Also: 0 &lt;= <paramref name="cscRowIndices"/>[i] &lt;
		/// <paramref name="order"/>.
		/// </param>
		/// <param name="cscColOffsets">
		/// Contains the index of the first entry of each column into the arrays <paramref name="cscValues"/> and
		/// <paramref name="cscRowIndices"/>. Its length must be <paramref name="order"/> + 1. The last entry must be
		/// <paramref name="numNonZerosUpper"/>.
		/// </param>
		/// <exception cref="IndefiniteMatrixException">Thrown if the original matrix is not positive definite.</exception>
		void Factorize(int order, int numNonZerosUpper, double[] cscValues, int[] cscRowIndices, int[] cscColOffsets);

		/// <summary>
		/// Performs the Cholesky factorization: A = L * L^T of a symmetric positive definite matrix A.
		/// Only the upper triangle of the original matrix is required and is provided in symmetric CSC format.
		/// </summary>
		/// <param name="matrix">The matrix in symmetric (only upper triangle) CSC format.</param>
		/// <exception cref="IndefiniteMatrixException">Thrown if the original matrix is not positive definite.</exception>
		void Factorize(SymmetricCscMatrix matrix);
	}
}
