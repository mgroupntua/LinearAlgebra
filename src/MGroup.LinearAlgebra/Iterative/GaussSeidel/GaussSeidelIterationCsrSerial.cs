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
		private CsrMatrix matrix;
		private int[] diagonalOffsets;
		private bool inactive = false;

		public GaussSeidelIterationCsrSerial()
		{
		}

		public void Dispose() 
		{
			this.diagonalOffsets = null;
			inactive = true;
		}

		public void GaussSeidelBackwardIteration(IVectorView rhsVector, IVector lhsVector)
		{
			CheckActive();

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
			CheckActive();

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

		public void Initialize(IMatrixView matrix) 
		{
			if (matrix is CsrMatrix csrMatrix)
			{
				Preconditions.CheckSquare(csrMatrix);
				this.matrix = csrMatrix;
				this.diagonalOffsets = SparseArrays.LocateDiagonalOffsets(
					matrix.NumRows, csrMatrix.RawRowOffsets, csrMatrix.RawColIndices);
				this.inactive = false;
			}
			else
			{
				throw new InvalidSparsityPatternException(this.GetType().Name + " can be used only for matrices in CSR format." +
					"Consider using the general " + typeof(GaussSeidelIterationGeneral).Name + " or another implementation.");
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
		private void CheckActive()
		{
			if (inactive)
			{
				throw new ObjectDisposedException(this.GetType().Name);
			}
		}

		public class Builder : IGaussSeidelIterationBuilder
		{
			public IGaussSeidelIteration Create() => new GaussSeidelIterationCsrSerial();
		}
	}
}
