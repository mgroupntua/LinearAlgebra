namespace MGroup.LinearAlgebra.Implementations.Managed
{
	using System;
	using System.Collections.Generic;
	using System.Text;

	using CSparse;
	using CSparse.Double;
	using CSparse.Ordering;

	using MGroup.LinearAlgebra.Reordering;

	/// <summary>
	/// Implementations of reordering algorithms provided by the CSparse.NET library.
	/// </summary>
	public class ManagedReorderingProvider : IReorderingProvider
	{
		/// <inheritdoc/>
		public (int[] permutation, ReorderingStatistics stats) AmdSymmetric(int order, int[] cscRowIndices, int[] cscColOffsets)
		{
			var stats = ReorderingStatistics.CreateUnknown(); // This implementation does not keep track

			var dummyCscValues = new double[cscRowIndices.Length]; //TODO: too expensive 
			var matrixCSparse = new SparseMatrix(order, order, dummyCscValues, cscRowIndices, cscColOffsets);
			int[] permutation = AMD.Generate<double>(matrixCSparse, ColumnOrdering.MinimumDegreeAtPlusA);

			// It is possible that CSparse.NET AMD algorithm returns more entries than the matrix order (so far I have found 1 
			// extra). In that case, make sure the first ones are valid and return only them.
			if (permutation.Length > order)
			{
				for (int i = order; i < permutation.Length; ++i)
				{
					if (permutation[i] < order)
					{
						throw new Exception(
							"Something went wrong during AMD. The permutation vector has more entries than the matrix order.");
					}
				}

				var permutationCorrected = new int[order];
				Array.Copy(permutation, permutationCorrected, order);
				return (permutationCorrected, stats);
			}
			else
			{
				return (permutation, stats);
			}
		}
	}
}
