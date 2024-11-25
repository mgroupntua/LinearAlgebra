////TODO: Creating a dummy variables array with as many entries as the non zero entries is expensive and not needed. I would be 
////      better off copying the AMD source code and avoiding that step.
////TODO: Find out what is going wrong when AMD returns a permutation with more entries than the matrix order.
//namespace MGroup.LinearAlgebra.Reordering
//{
//	using System;

//	using CSparse;
//	using CSparse.Double;
//	using CSparse.Ordering;

//	/// <summary>
//	/// Calculates a fill-reducing permutation for the rows/columns of a symmetric sparse matrix, using the Approximate Minimum 
//	/// Degree (AMD) ordering algorithm. The AMD implementation used is provided by the CSparse.NET library. For more 
//	/// information, see the AMD user guide, which is distributed as part of the SuiteSparse library.
//	/// </summary>
//	public class OrderingAmdCSparseNet : IReorderingAlgorithm
//	{
//		/// <inheritdoc/>
//		/// <remarks>The returned permutation is new-to-old.</remarks>
//		public (int[] permutation, bool oldToNew) FindPermutation(SparsityPatternSymmetric pattern)
//		{
//			int order = pattern.Order;
//			(int[] cscRowIndices, int[] cscColOffsets) = pattern.BuildSymmetricCSCArrays(true); //TODO: perhaps sorting is not needed here.

//			return FindPermutation(order, cscRowIndices, cscColOffsets);
//		}

//		/// <summary>
//		/// Finds a fill-reducting permutation for the rows/columns of a symmetric sparse matrix, described by its sparsity 
//		/// pattern. The returned permutation can be new-to-old or old-to-new.
//		/// </summary>
//		/// <param name="order">The number of rows (and columns) of the matrix.</param>
//		/// <param name="cscRowIndices">The CSC array that contains the row index of each non-zero entry.</param>
//		/// <param name="cscColOffsets">
//		/// The CSC array that contains the offset into <paramref name="cscRowIndices"/> of each column.
//		/// </param>
//		/// <returns>
//		/// permutation: An array containing the fill reducing permutation.
//		/// oldToNew: If false, then the permutation is defined as reordered[i] = original[permutation[i]].
//		///     If true, the permutation is defined as reordered[permutation[i]] = original[i].
//		/// </returns>
//		/// <remarks>The returned permutation is new-to-old.</remarks>
//		public (int[] permutation, bool oldToNew) FindPermutation(int order, int[] cscRowIndices, int[] cscColOffsets)
//		{
//			var dummyCscValues = new double[cscRowIndices.Length]; //TODO: too expensive 
//			var matrixCSparse = new SparseMatrix(order, order, dummyCscValues, cscRowIndices, cscColOffsets);
//			int[] permutation = AMD.Generate<double>(matrixCSparse, ColumnOrdering.MinimumDegreeAtPlusA);

//			// It is possible that CSparse.NET AMD algorithm returns more entries than the matrix order (so far I have found 1 
//			// extra). In that case, make sure the first ones are valid and return only them.
//			if (permutation.Length > order)
//			{
//				for (int i = order; i < permutation.Length; ++i)
//				{
//					if (permutation[i] < order)
//					{
//						throw new Exception(
//							"Something went wrong during AMD. The permutation vector has more entries than the matrix order.");
//					}
//				}

//				var permutationCorrected = new int[order];
//				Array.Copy(permutation, permutationCorrected, order);
//				return (permutationCorrected, false);
//			}
//			else
//			{
//				return (permutation, false);
//			}
//		}
//	}
//}
