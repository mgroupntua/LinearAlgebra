// TODO: also return the nonzeros after cholesky, flop count and other statistics
//TODO: the number of "dense" rows moved to the end is not reported by SuiteSparse and a dummy value (-1) is returned. Fix this.
namespace MGroup.LinearAlgebra.Implementations.NativeWin64.SuiteSparse
{
	using System;
	using System.Collections.Generic;
	using System.Text;

	using MGroup.LinearAlgebra.Exceptions;
	using MGroup.LinearAlgebra.Reordering;

	/// <summary>
	/// Implementations of reordering algorithms provided by the SuiteSparse library.
	/// </summary>
	public class SuiteSparseReorderingProvider : IReorderingProvider
	{
		/// <inheritdoc/>
		/// <exception cref="NativeLibException">Thrown if SuiteSparse dlls cannot be loaded or if AMD fails.</exception>
		public (int[] permutation, ReorderingStatistics stats) AmdSymmetric(int order, int[] cscRowIndices, int[] cscColOffsets)
		{
			var permutation = new int[order];
			IntPtr common = SuiteSparsePInvokes.CreateCommon(0, 0);
			if (common == IntPtr.Zero)
			{
				throw new NativeLibException("SuiteSparse (win64)", "Failed to initialize SuiteSparse.");
			}

			int numNonZerosUpper = cscRowIndices.Length;
			int status = SuiteSparsePInvokes.ReorderAMDUpper(order, numNonZerosUpper, cscRowIndices, cscColOffsets, permutation,
				out int nnzFactor, common);
			if (status == 0)
			{
				throw new NativeLibException("SuiteSparse (win64)", "AMD failed. This could be caused by the matrix being so " +
					"large it cannot be processed with the available memory.");
			}

			SuiteSparsePInvokes.DestroyCommon(ref common);
			return (permutation, new ReorderingStatistics(nnzFactor, -1));
		}

		/// <summary>
		/// Implements the Constrained Approximate Minimum Degree (AMD) fill-reducing ordering algorithm, by calling the 
		/// appropriate SuiteSparse library functions. CAMD is a variant of AMD that reorders the rows/columns of a matrix within
		/// each user defined group. On output, the rows/columns of each group are consecutive and each group is ordered 
		/// according to its user defined priority.
		/// The returned permutation is new-to-old, namely reordered[i] = original[permutation[i]].
		/// For more, see the CAMD user guide, which is distributed as part of the SuiteSparse library.
		/// </summary>
		/// <param name="order">The number of rows/columns of the symmetric matrix.</param>
		/// <param name="cscRowIndices">
		/// Row indices of the upper triangle entries of the symmetric matrix, in Compressed Sparse Columns format.
		/// All row indices of the same column must be sorted.
		/// </param>
		/// <param name="cscColOffsets">
		/// Column offsets of the upper triangle entries of the symmetric matrix, in Compressed Sparse Columns format.
		/// All column offsets must be sorted.
		/// </param>
		/// <param name="constraints">
		/// Array of length = order with ordering constraints.
		/// Its values must be 0 &lt;= <paramref name="constraints"/>[i] &lt; order.
		/// If <paramref name="constraints"/> = NULL, no constraints will be enforced.
		/// Example: <paramref name="constraints"/> = { 2, 0, 0, 0, 1 }. This means that indices 1, 2, 3 that have 
		/// <paramref name="constraints"/>[i] = 0, will be ordered before index 4 with <paramref name="constraints"/>[4] = 1,
		/// which will be ordered before index 0 with <paramref name="constraints"/>[0] = 2. Indeed for a certain pattern, 
		/// permutation = { 3, 2, 1, 4, 0 } (remember permutation is a new-to-old 
		/// mapping).
		/// </param>
		/// <param name="denseThreshold">
		/// A dense row/column in A + A^T can cause CAMD to spend significant time in ordering the matrix A.
		/// If <paramref name="denseThreshold"/> &gt;= 0, rows/columns with more than 
		/// <paramref name="denseThreshold"/> * sqrt(order(A)) entries are ignored during the ordering, and placed last in
		/// the output order. The default value of <paramref name="denseThreshold"/> is 10. If negative, no rows/columns are
		/// treated as dense. Rows/columns with 16 or fewer off-diagonal entries are never considered dense.
		/// WARNING: allowing dense rows/columns may violate the constraints (more tesing is needed).
		/// </param>
		/// <param name="aggressiveAbsorption">
		/// If true, aggressive absorption will be performed, which means that a prior element is absorbed into the current
		/// element if it is a subset of the current element, even if it is not adjacent to the current pivot element.
		/// This nearly always leads to a better ordering (because the approximate degrees are more accurate) and a lower
		/// execution time. However, there are cases where it can lead to a slightly worse ordering.
		/// </param>
		/// <returns>
		/// permutation: An array containing the new-to-old fill reducing permutation.
		/// stats: Measuments taken by SuiteSparse during the execution of AMD.
		/// </returns>
		/// <exception cref="ArgumentException">
		/// If <paramref name="order"/>, <paramref name="cscRowIndices"/> or <paramref name="cscColOffsets"/> do not describe 
		/// a valid symmetric matrix.
		/// </exception>
		/// <exception cref="NativeLibException">
		/// Thrown if SuiteSparse dlls cannot be loaded, or if there is not enough memory to allocate during CAMD.
		/// </exception>
		public (int[] permutation, ReorderingStatistics stats) Camd(int order, int[] cscRowIndices, int[] cscColOffsets,
			int[] constraints, int denseThreshold = -1, bool aggressiveAbsorption = true)
		{
			var permutation = new int[order];
			int status = SuiteSparsePInvokes.ReorderCAMD(order, cscRowIndices, cscColOffsets, constraints,
				denseThreshold, aggressiveAbsorption ? 1 : 0,
				permutation, out int nnzFactor, out int numMovedDense);

			// Error checking
			//if (status == 1)
			//{
			//	throw new SuiteSparseException("The matrix had unsorted columns, but was otherwise ok.");
			//}
			if (status == 2)
			{
				throw new ArgumentException("The arrays that describe the non zero pattern of the matrix were" +
					" invalid or the permutation buffer was NULL.");
			}
			else if (status == 3)
			{
				throw new NativeLibException("SuiteSparse (win64)", "Not enough memory could be allocated");
			}

			return (permutation, new ReorderingStatistics(nnzFactor, numMovedDense));
			//TODO: also return the nonzeros after cholesky, flop count and other statistics
		}
	}
}
