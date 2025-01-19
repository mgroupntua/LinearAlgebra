namespace MGroup.LinearAlgebra.Commons
{
	using System;
	using System.Collections.Generic;
	using System.Text;

	using MGroup.LinearAlgebra.Matrices;
	using MGroup.LinearAlgebra.Implementations;

	internal static class TransposeUtilities
	{
		internal static (TransposeMatrix transposeA, int lhsLength, int rhsLength) PrepareBlas(IIndexable2D matrix, bool transpose)
		{
			if (transpose)
			{
				return (TransposeMatrix.Transpose, matrix.NumRows, matrix.NumColumns);
			}
			else
			{
				return (TransposeMatrix.NoTranspose, matrix.NumColumns, matrix.NumRows);
			}
		}
	}
}
