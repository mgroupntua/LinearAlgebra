namespace MGroup.LinearAlgebra.SchurComplements
{
	using MGroup.LinearAlgebra.Matrices;
	using MGroup.LinearAlgebra.Triangulation;
	using MGroup.LinearAlgebra.Vectors;

	/// <summary>
	/// Calculates the Schur complement, where the original matrix A is partitioned into 4 submatrices A = [A00, A01; A10, A11]
	/// (in Matlabb notation). Specifically A is symmetric, A00 is in symmetric, packed, upper, column-major format
	/// (<see cref="SymmetricMatrix"/>), A01 is in CSR format (<see cref="CsrMatrix"/>), A10 is implied as the transpose of A01
	/// and A11 is provided through some factorization.
	/// </summary>
	public static class SchurComplementPckCsrSymCsc
	{
		/// <summary>
		/// Calculates the Schur complement of A/A11 = S = A00 - A01 * inv(A11) * A01^T, where A = [A00, A01; A10^T, A11]
		/// (in Matlab notation).
		/// </summary>
		/// <param name="A01">It must be in typical CSR format, where the rows must contain columns in ascending order.</param>
		/// <remarks>
		/// This method constructs inv(A11) * A10 one column at a time and uses that column to calculate the superdiagonal
		/// entries of the corresponding column of A01 * inv(A11) * A10.
		/// </remarks>
		public static SymmetricMatrix CalcSchurComplement(SymmetricMatrix A00, CsrMatrix A01, ITriangulation inverseA11)
		{
			var S = SymmetricMatrix.CreateZero(A00.Order);
			CalcSchurComplement(A00, A01, inverseA11, S);
			return S;
		}

		/// <summary>
		/// Calculates the Schur complement of A/A11 = S = A00 - A01 * inv(A11) * A01^T, where A = [A00, A01; A10^T, A11]
		/// (in Matlab notation).
		/// </summary>
		/// <param name="A01">It must be in typical CSR format, where the rows must contain columns in ascending order.</param>
		/// <remarks>
		/// This method constructs inv(A11) * A10 one column at a time and uses that column to calculate the superdiagonal
		/// entries of the corresponding column of A01 * inv(A11) * A10.
		/// </remarks>
		public static void CalcSchurComplement(SymmetricMatrix A00, CsrMatrix A01, ITriangulation inverseA11, 
			SymmetricMatrix result)
		{ 
			//TODO: Unfortunately this cannot take advantage of MKL for CSR * matrix or even CSR * vector. Instead there should
			////	be faster implementations of this whole class 
			for (int j = 0; j < A01.NumRows; ++j)
			{
				// Column j of A10 = row j of A01
				Vector colA10 = A01.GetRow(j);

				// column j of (inv(A11) * A10) = inv(A00) * column j of A10
				double[] colInvA11TimesA10 = inverseA11.SolveLinearSystem(colA10).RawData;

				// column j of (A01 * inv(11) * A10) = A01 * column j of (inv(A00) * A10)
				// However we only need the superdiagonal part of this column. 
				// Thus we only multiply the rows i of A01 with i <= j. 
				for (int i = 0; i <= j; ++i)
				{
					// Perform the subtraction S = A00 - (A10^T * inv(A11) * A10) for the current (i, j)
					double dot = MultiplyRowTimesVector(A01, i, colInvA11TimesA10);
					int indexS = i + (j * (j + 1)) / 2;
					result.RawValues[indexS] = A00.RawValues[indexS] - dot;
				}
			}
		}

		private static double MultiplyRowTimesVector(CsrMatrix csr, int rowIdx, double[] vector)
		{
			double[] values = csr.RawValues;
			int[] colIndices = csr.RawColIndices;
			int[] rowOffsets = csr.RawRowOffsets;
			double dot = 0.0;
			int start = rowOffsets[rowIdx]; //inclusive
			int end = rowOffsets[rowIdx + 1]; //exclusive
			for (int k = start; k < end; ++k)
			{
				dot += values[k] * vector[colIndices[k]];
			}

			return dot;
		}
	}
}
