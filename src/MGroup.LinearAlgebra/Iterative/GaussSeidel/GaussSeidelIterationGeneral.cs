namespace MGroup.LinearAlgebra.Iterative.GaussSeidel
{
	using System;
	using System.Collections.Generic;
	using System.Text;

	using MGroup.LinearAlgebra.Commons;
	using MGroup.LinearAlgebra.Matrices;
	using MGroup.LinearAlgebra.Vectors;

	public class GaussSeidelIterationGeneral : IGaussSeidelIteration
	{
		private readonly IMatrixView matrix;

		public GaussSeidelIterationGeneral(IMatrixView matrix)
		{
			Preconditions.CheckSquare(matrix);
			this.matrix = matrix;
			this.SystemSize = matrix.NumRows;
		}

		public int SystemSize { get; }

		public void Dispose() { }

		public void GaussSeidelBackwardIteration(IVectorView rhsVector, IVector lhsVector)
		{
			IVector x = lhsVector;
			IVectorView b = rhsVector;
			Preconditions.CheckSquareLinearSystemDimensions(matrix, x, b);
			int n = matrix.NumRows;
			for (int i = n - 1; i >= 0; --i)
			{
				double sum = b[i];
				int j;
				for (j = n - 1; j > i; --j)
				{
					sum -= matrix[i, j] * x[j];
				}

				double diagEntry = matrix[i, i];

				for (j = i - 1; j >= 0; --j)
				{
					sum -= matrix[i, j] * x[j];
				}

				x.Set(i, sum / diagEntry);
			}
		}

		public void GaussSeidelForwardIteration(IVectorView rhsVector, IVector lhsVector)
		{
			IVector x = lhsVector;
			IVectorView b = rhsVector;
			Preconditions.CheckSquareLinearSystemDimensions(matrix, x, b);
			int n = matrix.NumRows;
			for (int i = 0; i < n; ++i)
			{
				double sum = b[i];
				int j;
				for (j = 0; j < i; ++j)
				{
					sum -= matrix[i, j] * x[j];
				}

				double diagEntry = matrix[i, i];

				for (j = i + 1; j < n; ++j)
				{
					sum -= matrix[i, j] * x[j];
				}

				x.Set(i, sum / diagEntry);
			}
		}

		public void Initialize() { }
	}
}
