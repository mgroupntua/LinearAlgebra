namespace MGroup.LinearAlgebra.Reordering
{
	using System;
	using System.Collections.Generic;
	using System.Text;

	using MGroup.LinearAlgebra.Implementations;
	using MGroup.LinearAlgebra.Implementations.Managed;
	using MGroup.LinearAlgebra.Matrices.Builders;

	using static MGroup.LinearAlgebra.LibrarySettings;

	/// <summary>
	/// Calculates a fill-reducing permutation for the rows/columns of a symmetric sparse matrix, using the Approximate Minimum
	/// Degree (AMD) ordering algorithm.
	/// For more information, see the AMD user guide, which is distributed as part of the SuiteSparse library.
	/// </summary>
	public class AmdSymmetricOrdering : IReorderingAlgorithm
	{

		public AmdSymmetricOrdering()
		{
		}

		/// <inheritdoc/>
		/// <remarks>The returned permutation is always new-to-old when using AMD.</remarks>
		public (int[] permutation, bool oldToNew) FindPermutation(SparsityPatternSymmetric pattern)
		{
			(int[] rowIndices, int[] colOffsets) = pattern.BuildSymmetricCSCArrays(sortRowsOfEachCol: true);
			(int[] permutation, _) = GlobalProvider.Reordering.AmdSymmetric(pattern.Order, rowIndices, colOffsets);
			return (permutation, false);
		}

		/// <inheritdoc/>
		/// <remarks>The returned permutation is always new-to-old when using AMD.</remarks>
		public (int[] permutation, bool oldToNew) FindPermutation(int order, int[] cscRowIndices, int[] cscColOffsets)
		{
			(int[] permutation, _) = GlobalProvider.Reordering.AmdSymmetric(order, cscRowIndices, cscColOffsets);
			return (permutation, false);
		}

		/// <summary>
		/// Finds a fill-reducting permutation for the rows/columns of a symmetric sparse matrix in DOK format.
		/// </summary>
		/// <remarks>The returned permutation is always new-to-old when using AMD.</remarks>
		/// <param name="dok">The symmetric sparse matrix in DOK format.</param>
		/// <returns>
		/// permutation: An array containing the fill reducing permutation.
		/// stats: Measuments taken during the execution of the reordering algorithm.
		/// </returns>
		public (int[] permutation, bool oldToNew, ReorderingStatistics stats) FindPermutation(DokSymmetric dok)
		{
			(double[] values, int[] rowIndices, int[] colOffsets) = dok.BuildSymmetricCscArrays(true);
			(int[] permutation, ReorderingStatistics stats) =
				GlobalProvider.Reordering.AmdSymmetric(dok.NumColumns, rowIndices, colOffsets);
			return (permutation, false, stats);
		}
	}
}
