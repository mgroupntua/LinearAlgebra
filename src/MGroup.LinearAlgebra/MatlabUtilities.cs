using System;
using System.Collections.Generic;
using System.Text;
using MGroup.LinearAlgebra.Matrices;
using MGroup.LinearAlgebra.Vectors;

namespace MGroup.LinearAlgebra
{
	public static class MatlabUtilities
	{
		public static int[] Find(this Vector vector)
		{
			var result = new List<int>(vector.Length);
			for (int i = 0; i < vector.Length; ++i)
			{
				if (vector[i] != 0)
				{
					result.Add(i);
				}
			}
			return result.ToArray();
		}

		public static int[] Find(this Vector vector, Predicate<double> predicate)
		{
			var result = new List<int>(vector.Length);
			for (int i = 0; i < vector.Length; ++i)
			{
				if (predicate(vector[i]))
				{
					result.Add(i);
				}
			}
			return result.ToArray();
		}

		/// <summary>
		/// Creates and returns an m-by-n sparse matrix S by taking the columns of this matrix B and placing them along the 
		/// diagonals specified by d
		/// </summary>
		/// <param name="B"></param>
		/// <param name="d"></param>
		/// <param name="m"></param>
		/// <param name="n"></param>
		/// <returns></returns>
		public static Matrix Spdiags(this Matrix B, int[] d, int m, int n)
		{
			if (B.NumColumns != d.Length)
			{
				throw new ArgumentException();
			}

			var S = Matrix.CreateZero(m, n);
			if (m >= n)
			{
				int numActiveRows = Math.Min(B.NumRows, n);
				for (int j = 0; j < d.Length; ++j)
				{
					int rowOffset = -d[j];
					int iStart, iEnd;
					if (d[j] >= 0)
					{
						iStart = d[j];
						iEnd = numActiveRows - 1;
					}
					else
					{
						iStart = 0;
						iEnd = m - 1 - rowOffset;
					}
					for (int i = iStart; i <= iEnd; ++i)
					{
						int row = rowOffset + i;
						int col = i;
						S[row, col] = B[i, j];
					}
				}
			}
			else
			{
				int numActiveRows = Math.Min(B.NumRows, m);
				for (int j = 0; j < d.Length; ++j)
				{
					int colOffset = +d[j];
					int iStart, iEnd;
					if (d[j] > 0)
					{
						iStart = 0;
						iEnd = n - 1 - colOffset;
					}
					else
					{
						iStart = -d[j];
						iEnd = numActiveRows - 1;
					}
					for (int i = iStart; i <= iEnd; ++i)
					{
						int row = i;
						int col = colOffset + i;
						S[row, col] = B[i, j];
					}
				}
			}

			return S;
		}
	}
}
