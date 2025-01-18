namespace MGroup.LinearAlgebra.Triangulation
{
	using System;
	using System.Collections.Generic;
	using System.Text;

	using MGroup.LinearAlgebra.Exceptions;
	using MGroup.LinearAlgebra.Matrices;

	public interface ILUCscFactorization : ITriangulation
	{
		/// <summary>
		/// Performs the LU factorization: A = L * U of a ssquare matrix A.  The matrix A is provided in CSC format described by
		/// <paramref name="cscValues"/>, <paramref name="cscRowIndices"/> and <paramref name="cscColOffsets"/>.
		/// </summary>
		/// <param name="order">The number of rows/columns of the square matrix.</param>
		/// <param name="numNonZeros">The number of explicitly stored entries of the matrix before the factorization.</param>
		/// <param name="cscValues">
		/// Contains the non-zero entries. Its length must be equal to <paramref name="numNonZeros"/>.
		/// The non-zero entries of each row must appear consecutively in <paramref name="cscValues"/>. They should also be
		/// sorted in increasing order of their row indices, to speed up subsequent the factorization.
		/// </param>
		/// <param name="cscRowIndices">
		/// Contains the row indices of the non-zero entries. Its length must be equal to <paramref name="numNonZeros"/>.
		/// There is an 1 to 1 matching between these two arrays: <paramref name="cscRowIndices"/>[i] is the row index of the
		/// entry <paramref name="cscValues"/>[i]. Also: 0 &lt;= <paramref name="cscRowIndices"/>[i] &lt;
		/// <paramref name="order"/>.
		/// </param>
		/// <param name="cscColOffsets">
		/// Contains the index of the first entry of each column into the arrays <paramref name="cscValues"/> and
		/// <paramref name="cscRowIndices"/>. Its length must be <paramref name="order"/> + 1. The last entry must be
		/// <paramref name="numNonZeros"/>.
		/// </param>
		/// <param name="pivotTolerance">The partial pivoting tolerance (from 0.0 to 1.0).</param>
		/// <exception cref="SingularMatrixException">Thrown if the original matrix is not invertible.</exception>
		public void Factorize(
			int order, int numNonZeros, double[] cscValues, int[] cscRowIndices, int[] cscColOffsets, double pivotTolerance);

		/// <summary>
		/// Performs the LU factorization: A = L * U of a square matrix A. The matrix A is provided in CSC format.
		/// </summary>
		/// <param name="matrix">The matrix in CSC format. Must be a square matrix.</param>
		/// <param name="pivotTolerance">The partial pivoting tolerance (from 0.0 to 1.0).</param>
		/// <exception cref="SingularMatrixException">Thrown if the original matrix is not invertible.</exception>
		public void Factorize(CscMatrix matrix, double pivotTolerance);
	}
}
