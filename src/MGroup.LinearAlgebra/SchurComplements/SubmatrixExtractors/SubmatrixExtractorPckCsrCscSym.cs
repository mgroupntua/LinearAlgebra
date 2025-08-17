namespace MGroup.LinearAlgebra.SchurComplements.SubmatrixExtractors
{
	using System;

	using MGroup.LinearAlgebra.Commons;
	using MGroup.LinearAlgebra.Matrices;
	using MGroup.LinearAlgebra.SchurComplements.IntegerMatrices;

	/// <summary>
	/// <inheritdoc path="/summary"/>
	/// Submatrix A00 is in packed, upper, column-major format.
	/// </summary>
	public class SubmatrixExtractorPckCsrCscSym : SubmatrixExtractorCscSymBase
	{
		/// <summary>
		/// A00 of A = [A00, A01; A10, A11] (using Matlab notation).
		/// </summary>
		public SymmetricMatrix Submatrix00 { get; private set; }

		/// <inheritdoc/>
		public override void Clear()
		{
			base.Clear();
			Submatrix00 = null;
		}

		/// <summary>
		/// Calculates the submatrices of A=<paramref name="originalMatrix"/>, such that A = [A00, A01; A10, A11]
		/// (in Matlab notation). In this partition, the rows and columns that correspond to groups 0 and 1 are defined in
		/// <paramref name="indicesGroup0"/> and <paramref name="indicesGroup1"/> respectively. The resulting submatrices are
		/// stored inside the properties of this object.
		/// </summary>
		/// <param name="originalMatrix">The original matrix A that will be partitioned.</param>
		/// <param name="indicesGroup0">
		/// The rows and columns of A that will end up as rows and columns of A00.
		/// They do not need to be contiguous or in any order.
		/// </param>
		/// <param name="indicesGroup1">
		/// The rows and columns of A that will end up as rows and columns of A11.
		/// They do not need to be contiguous or in any order.
		/// </param>
		/// <exception cref="InvalidOperationException">
		/// Thrown if a different decomposition was performed before (matrix with different sparsity pattern or different
		/// index groups), without clearing the stored data of this object.
		/// </exception>
		public void ExtractSubmatrices(SymmetricCscMatrix originalMatrix, int[] indicesGroup0, int[] indicesGroup1)
		{
			if (this.originalMatrix == null)
			{
				ExtractMaps(originalMatrix, indicesGroup0, indicesGroup1);
				this.originalMatrix = originalMatrix;
				this.indicesGroup0 = indicesGroup0;
				this.indicesGroup1 = indicesGroup1;
			}
			else
			{
				if ((this.originalMatrix != originalMatrix)
					|| (this.indicesGroup0 != indicesGroup0)
					|| (this.indicesGroup1 != indicesGroup1))
				{
					throw new InvalidOperationException("The submatrix pattern stored corresponds to a different original matrix"
						+ " or to different subsets of rows/columns. Clear this object first");
				}
			}

			mapper00.CopyValuesArrayToSubmatrix(originalMatrix.RawValues, Submatrix00.RawValues);
			mapper01.CopyValuesArrayToSubmatrix(originalMatrix.RawValues, Submatrix01.RawValues);
			mapper11.CopyValuesArrayToSubmatrix(originalMatrix.RawValues, Submatrix11.RawValues);
		}

		private void ExtractMaps(SymmetricCscMatrix originalMatrix, int[] indicesGroup0, int[] indicesGroup1)
		{
			int n0 = indicesGroup0.Length;
			int n1 = indicesGroup1.Length;
			var submatrix00 = IntSymMatrixColMajor.CreateZero(n0);
			ArrayUtilities.MemSet(submatrix00.RawValues, -1); // Later, -1 index will be overwritten by all entries, except explicit zeros of A00 
			var submatrix01 = IntDokRowMajor.CreateZero(n0, n1);
			var submatrix11 = IntDokSymColMajor.CreateZero(n1);

			// Original matrix indices to submatrix indices. Group 0 indices i0 are stored as: -i0-1.
			// Group 1 indices are stored as they are
			int[] originalToSubIndices = MapOriginalToSubmatrixIndices(originalMatrix.NumRows, indicesGroup0, indicesGroup1);

			// Iterate the non zero values array of the original sparse matrix
			for (int j = 0; j < originalMatrix.NumColumns; ++j)
			{
				int start = originalMatrix.RawColOffsets[j];
				int end = originalMatrix.RawColOffsets[j + 1];
				int subJ = originalToSubIndices[j];
				if (subJ < 0) // j belongs to group 0 indices
				{
					subJ = -subJ - 1; // Mapping to [0, oo) for group 0 indices
					for (int t = start; t < end; ++t)
					{
						int i = originalMatrix.RawRowIndices[t];
						int subI = originalToSubIndices[i];
						if (subI < 0) // (i, j) belongs to submatrix A00
						{
							subI = -subI- 1; // Mapping to [0, oo) for group 0 indices
							submatrix00[subI, subJ] = t;
							submatrix00[subJ, subI] = t;
						}
						else // (i, j) belongs to submatrix A10 or equivalently (j, i) belongs to submatrix A01
						{
							submatrix01[subJ, subI] = t;
						}
					}
				}
				else // j belongs to group 1 indices
				{
					for (int t = start; t < end; ++t)
					{
						int i = originalMatrix.RawRowIndices[t];
						int subI = originalToSubIndices[i];
						if (subI < 0) // (i, j) belongs to submatrix A01
						{
							subI = -subI- 1; // Mapping to [0, oo) for group 0 indices
							submatrix01[subI, subJ] = t;
						}
						else // (i, j) belongs to submatrix A11
						{
							submatrix11[subI, subJ] = t;
						}
					}
				}
			}

			// Finalize the data structures required to represent the submatrices
			// A00 packed, col major. Zero entries that were not stored in A will have a -1 index into A.RawValues.
			this.Submatrix00 = SymmetricMatrix.CreateZero(n0);
			this.mapper00 = new SparseToDenseValuesArrayMapper(submatrix00.RawValues);

			// A01 CSR
			(int[] submatrixToOriginalValues01, int[] colIndices01, int[] rowOffsets01) = submatrix01.BuildCsrArrays();
			this.Submatrix01 = CsrMatrix.CreateFromArrays(
				n0, n1, new double[colIndices01.Length], colIndices01, rowOffsets01, false);
			this.mapper01 = new SameSparsityValuesArrayMapper(submatrixToOriginalValues01);

			// A11 CSC upper triangle
			(int[] submatrixToOriginalValues11, int[] rowIndices11, int[] colOffsets11) = submatrix11.BuildCscArrays();
			this.Submatrix11 = SymmetricCscMatrix.CreateFromArrays(
				n1, new double[rowIndices11.Length], rowIndices11, colOffsets11, false);
			this.mapper11 = new SameSparsityValuesArrayMapper(submatrixToOriginalValues11);
		}
	}
}
