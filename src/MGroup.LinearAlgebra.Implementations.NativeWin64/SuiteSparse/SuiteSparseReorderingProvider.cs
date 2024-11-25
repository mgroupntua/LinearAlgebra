// TODO: also return the nonzeros after cholesky, flop count and other statistics
//TODO: the number of "dense" rows moved to the end is not reported by SuiteSparse and a dummy value (-1) is returned. Fix this.
namespace MGroup.LinearAlgebra.Implementations.NativeWin64
{
	using System;
	using System.Collections.Generic;
	using System.Text;

	using MGroup.LinearAlgebra.Exceptions;
	using MGroup.LinearAlgebra.Implementations.PInvoke;
	using MGroup.LinearAlgebra.Reordering;

	/// <summary>
	/// Implementations of reordering algorithms provided by the SuiteSparse library.
	/// </summary>
	public class SuiteSparseReorderingProvider : IReorderingProvider
	{
		/// <inheritdoc/>
		/// <exception cref="SuiteSparseException">Thrown if SuiteSparse dlls cannot be loaded or if AMD fails.</exception>
		public (int[] permutation, ReorderingStatistics stats) AmdSymmetric(int order, int[] cscRowIndices, int[] cscColOffsets)
		{
			var permutation = new int[order];
			IntPtr common = SuiteSparsePInvokes.CreateCommon(0, 0);
			if (common == IntPtr.Zero)
			{
				throw new SuiteSparseException("Failed to initialize SuiteSparse.");
			}

			int numNonZerosUpper = cscRowIndices.Length;
			int status = SuiteSparsePInvokes.ReorderAMDUpper(order, numNonZerosUpper, cscRowIndices, cscColOffsets, permutation,
				out int nnzFactor, common);
			if (status == 0)
			{
				throw new SuiteSparseException("AMD failed. This could be caused by the matrix being so large it"
					+ " cannot be processed with the available memory.");
			}

			SuiteSparsePInvokes.DestroyCommon(ref common);
			return (permutation, new ReorderingStatistics(nnzFactor, -1));
		}
	}
}
