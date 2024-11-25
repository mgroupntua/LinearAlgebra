namespace MGroup.LinearAlgebra.Reordering
{
	/// <summary>
	/// Calculates fill-reducting permutations for the rows/columns of a symmetric sparse matrix.
	/// </summary>
	public interface IReorderingAlgorithm
	{
		/// <summary>
		/// Finds a fill-reducting permutation for the rows/columns of a symmetric sparse matrix, described by its sparsity 
		/// pattern. The returned permutation can be new-to-old or old-to-new.
		/// </summary>
		/// <param name="pattern">The indices of the non-zero entries of a symmetric matrix.</param>
		/// <returns>
		/// permutation: An array containing the fill reducing permutation.
		/// oldToNew: If false, then the permutation is defined as reordered[i] = original[permutation[i]].
		///     If true, the permutation is defined as reordered[permutation[i]] = original[i].
		/// </returns>
		(int[] permutation, bool oldToNew) FindPermutation(SparsityPatternSymmetric pattern);

		/// <summary>
		/// Finds a fill-reducting permutation for the rows/columns of a symmetric sparse matrix, described by its sparsity
		/// pattern. The returned permutation can be new-to-old or old-to-new.
		/// </summary>
		/// <param name="order">The number of rows (and columns) of the matrix.</param>
		/// <param name="cscRowIndices">
		/// Array that contains the row indices of each non-zero entry of the upper triangle, in symmetric CSC format.
		/// </param>
		/// <param name="cscColOffsets">
		/// Array that contains the offsets into <paramref name="cscRowIndices"/> of each column, in symmetric CSC format.
		/// </param>
		/// <returns>
		/// permutation: An array containing the fill reducing permutation.
		/// oldToNew: If false, then the permutation is defined as reordered[i] = original[permutation[i]].
		///     If true, the permutation is defined as reordered[permutation[i]] = original[i].
		/// </returns>
		(int[] permutation, bool oldToNew) FindPermutation(int order, int[] cscRowIndices, int[] cscColOffsets);
	}
}
