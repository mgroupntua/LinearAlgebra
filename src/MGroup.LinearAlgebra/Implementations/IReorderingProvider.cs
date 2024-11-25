namespace MGroup.LinearAlgebra.Implementations
{
	using System;
	using System.Collections.Generic;
	using System.Text;

	using MGroup.LinearAlgebra.Reordering;

	public interface IReorderingProvider
	{
		/// <summary>
		/// Apply the Approximate Minimum Degree (AMD) algorithm to find a fill-reducing permutation of the rows and columns.
		/// AMD's permutation is always new-to-old, namely reordered[i] = original[permutation[i]].
		/// </summary>
		/// <param name="order">The number of rows/columns of the matrix.</param>
		/// <param name="cscRowIndices">The row indices in symmetric CSC format of the upper triangle of the matrix.</param>
		/// <param name="cscColOffsets">The column offsets in symmetric CSC format of the upper triangle of the matrix.</param>
		/// <returns>permutation: a new-to-old mapping, stats: statistics about the execution of the algorithm.</returns>
		(int[] permutation, ReorderingStatistics stats) AmdSymmetric(int order, int[] cscRowIndices, int[] cscColOffsets);
	}
}
