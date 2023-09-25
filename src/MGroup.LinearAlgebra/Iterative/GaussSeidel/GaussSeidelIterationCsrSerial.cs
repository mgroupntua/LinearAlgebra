using System;
using System.Collections.Generic;
using System.Runtime.CompilerServices;
using System.Text;

using MGroup.LinearAlgebra.Commons;
using MGroup.LinearAlgebra.Exceptions;
using MGroup.LinearAlgebra.Matrices;
using MGroup.LinearAlgebra.Vectors;

namespace MGroup.LinearAlgebra.Iterative.GaussSeidel
{
	public class GaussSeidelIterationCsrSerial : IGaussSeidelIteration
	{
		private readonly CsrMatrix matrix;
		private int[] diagonalOffsets;
		private bool disposed = false;

		public int SystemSize { get; }

		public GaussSeidelIterationCsrSerial(CsrMatrix matrix)
		{
			Preconditions.CheckSquare(matrix);
			this.matrix = matrix;
			this.SystemSize = matrix.NumRows;
		}

		public void Dispose() 
		{
			this.diagonalOffsets = null;
			disposed = true;
		}

		public void GaussSeidelBackwardIteration(IVectorView rhsVector, IVector lhsVector)
		{
			CheckDisposed();

			if ((lhsVector is Vector lhsDense) && (rhsVector is Vector rhsDense))
			{
				BackwardIteration(matrix.NumRows, matrix.RawValues, matrix.RawRowOffsets, matrix.RawColIndices, diagonalOffsets,
					rhsDense.RawData, lhsDense.RawData);
			}
			else
			{
				double[] lhs = lhsVector.CopyToArray();
				double[] rhs = rhsVector.CopyToArray();
				BackwardIteration(matrix.NumRows, matrix.RawValues, matrix.RawRowOffsets, matrix.RawColIndices, diagonalOffsets, 
					rhs, lhs);
				lhsVector.CopyFrom(Vector.CreateFromArray(lhs));
			}
		}

		public void GaussSeidelForwardIteration(IVectorView rhsVector, IVector lhsVector)
		{
			CheckDisposed();

			if ((lhsVector is Vector lhsDense) && (rhsVector is Vector rhsDense))
			{
				ForwardIteration(matrix.NumRows, matrix.RawValues, matrix.RawRowOffsets, matrix.RawColIndices, diagonalOffsets,
					rhsDense.RawData, lhsDense.RawData);
			}
			else
			{
				double[] lhs = lhsVector.CopyToArray();
				double[] rhs = rhsVector.CopyToArray();
				ForwardIteration(matrix.NumRows, matrix.RawValues, matrix.RawRowOffsets, matrix.RawColIndices, diagonalOffsets, 
					rhs, lhs);
				lhsVector.CopyFrom(Vector.CreateFromArray(lhs));
			}
		}

		public void Initialize() 
		{
			CheckDisposed();

			int n = matrix.NumRows;
			int[] rowOffets = matrix.RawRowOffsets;
			int[] colIndices = matrix.RawColIndices;
			diagonalOffsets = new int[n];
			for (int i = 0; i < n; ++i)
			{
				int rowStart = rowOffets[i]; // inclusive
				int rowEnd = rowOffets[i + 1]; // exclusive

				//TODO: optimizations: bisection, start from the end of the row if row > n/2, etc.
				bool isDiagonalZero = true;
				for (int k = rowStart; k < rowEnd; ++k) 
				{ 
					int j = colIndices[k];
					if (j == i)
					{
						diagonalOffsets[i] = k;
						isDiagonalZero = false;
						break;
					}
				}

				if (isDiagonalZero)
				{
					throw new InvalidSparsityPatternException(
						$"Found diagonal entry at ({i},{i}). Gauss Seidel cannot be performed");
				}
			}
		}

		private static void BackwardIteration(int matrixOrder, double[] csrValues, int[] csrRowOffsets, int[] csrColIndices,
			int[] diagOffsets, double[] rhs, double[] lhs)
		{
			// Do not read the vector and each matrix row backward. It destroys caching. And in any case, we do csr_row * vector,
			// thus the order of operations for the dot product do not matter. What does matter is starting from the last row. 
			int n = matrixOrder;
			for (int i = n - 1; i >= 0; --i)
			{
				double sum = rhs[i];
				int rowStart = csrRowOffsets[i]; // inclusive
				int rowEnd = csrRowOffsets[i + 1]; // exclusive
				int diagOffset = diagOffsets[i];

				for (int k = rowStart; k < diagOffset; ++k)
				{
					sum -= csrValues[k] * lhs[csrColIndices[k]];
				}

				double diagEntry = csrValues[diagOffset];

				for (int k = diagOffset + 1; k < rowEnd; ++k)
				{
					sum -= csrValues[k] * lhs[csrColIndices[k]];
				}
				lhs[i] = sum / diagEntry;
			}
		}

		private static void BackwardIterationBasic(int matrixOrder, double[] csrValues, int[] csrRowOffsets,
			int[] csrColIndices, double[] rhs, double[] lhs)
		{
			int n = matrixOrder;
			for (int i = n - 1; i >= 0; --i)
			{
				double sum = rhs[i];
				double diagEntry = 0;

				int rowStart = csrRowOffsets[i]; // inclusive
				int rowEnd = csrRowOffsets[i + 1]; // exclusive
				for (int k = rowEnd - 1; k >= rowStart; --k)
				{
					int j = csrColIndices[k];
					if (j == i)
					{
						diagEntry = csrValues[k];
					}
					else
					{
						sum -= csrValues[k] * lhs[j];
					}
				}
				lhs[i] = sum / diagEntry;
			}
		}

		private static void ForwardIteration(int matrixOrder, double[] csrValues, int[] csrRowOffsets, int[] csrColIndices,
			int[] diagOffsets, double[] rhs, double[] lhs)
		{
			int n = matrixOrder;
			for (int i = 0; i < n; ++i)
			{
				double sum = rhs[i];
				int rowStart = csrRowOffsets[i]; // inclusive
				int rowEnd = csrRowOffsets[i + 1]; // exclusive
				int diagOffset = diagOffsets[i];
				
				for (int k = rowStart; k < diagOffset; ++k)
				{
					sum -= csrValues[k] * lhs[csrColIndices[k]];
				}

				double diagEntry = csrValues[diagOffset];

				for (int k = diagOffset + 1; k < rowEnd; ++k)
				{
					sum -= csrValues[k] * lhs[csrColIndices[k]];
				}
				lhs[i] = sum / diagEntry;
			}
		}

		private static void ForwardIterationBasic(int matrixOrder, double[] csrValues, int[] csrRowOffsets, int[] csrColIndices,
			double[] rhs, double[] lhs)
		{
			int n = matrixOrder;
			for (int i = 0; i < n; ++i)
			{
				double sum = rhs[i];
				double diagEntry = 0;

				int rowStart = csrRowOffsets[i]; // inclusive
				int rowEnd = csrRowOffsets[i + 1]; // exclusive
				for (int k = rowStart; k < rowEnd; ++k)
				{
					int j = csrColIndices[k];
					if (j == i)
					{
						diagEntry = csrValues[k];
					}
					else
					{
						sum -= csrValues[k] * lhs[j];
					}
				}
				lhs[i] = sum / diagEntry;
			}
		}

		[MethodImpl(MethodImplOptions.AggressiveInlining)]
		private void CheckDisposed()
		{
			if (disposed)
			{
				throw new ObjectDisposedException(this.GetType().Name);
			}
		}
	}
}
