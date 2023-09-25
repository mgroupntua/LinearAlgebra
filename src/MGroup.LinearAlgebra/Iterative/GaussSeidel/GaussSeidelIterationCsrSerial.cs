namespace MGroup.LinearAlgebra.Iterative.GaussSeidel
{
	using System;
	using System.Collections.Generic;
	using System.Text;

	using MGroup.LinearAlgebra.Commons;
	using MGroup.LinearAlgebra.Matrices;
	using MGroup.LinearAlgebra.Vectors;

	public class GaussSeidelIterationCsrSerial : IGaussSeidelIteration
	{
		private readonly CsrMatrix matrix;

		public int SystemSize { get; }

		public GaussSeidelIterationCsrSerial(CsrMatrix matrix)
		{
			Preconditions.CheckSquare(matrix);
			this.matrix = matrix;
			this.SystemSize = matrix.NumRows;
		}

		public void Dispose() 
		{ 
		
		}

		public void GaussSeidelBackwardIteration(IVectorView rhsVector, IVector lhsVector)
		{
			if ((lhsVector is Vector lhsDense) && (rhsVector is Vector rhsDense))
			{
				BackwardIteration(matrix.NumRows, matrix.RawValues, matrix.RawRowOffsets, matrix.RawColIndices, 
					rhsDense.RawData, lhsDense.RawData);
			}
			else
			{
				double[] lhs = lhsVector.CopyToArray();
				double[] rhs = rhsVector.CopyToArray();
				BackwardIteration(matrix.NumRows, matrix.RawValues, matrix.RawRowOffsets, matrix.RawColIndices, rhs, lhs);
				lhsVector.CopyFrom(Vector.CreateFromArray(lhs));
			}
		}

		public void GaussSeidelForwardIteration(IVectorView rhsVector, IVector lhsVector)
		{
			if ((lhsVector is Vector lhsDense) && (rhsVector is Vector rhsDense))
			{
				ForwardIteration(matrix.NumRows, matrix.RawValues, matrix.RawRowOffsets, matrix.RawColIndices,
					rhsDense.RawData, lhsDense.RawData);
			}
			else
			{
				double[] lhs = lhsVector.CopyToArray();
				double[] rhs = rhsVector.CopyToArray();
				ForwardIteration(matrix.NumRows, matrix.RawValues, matrix.RawRowOffsets, matrix.RawColIndices, rhs, lhs);
				lhsVector.CopyFrom(Vector.CreateFromArray(lhs));
			}
		}

		public void Initialize() 
		{ 

		}

		private static void BackwardIteration(int matrixOrder, double[] csrValues, int[] csrRowOffsets,
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
	}
}
