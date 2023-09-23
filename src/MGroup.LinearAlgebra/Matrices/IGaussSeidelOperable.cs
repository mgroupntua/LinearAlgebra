namespace MGroup.LinearAlgebra.Matrices
{
	using System;
	using System.Collections.Generic;
	using System.Text;

	using MGroup.LinearAlgebra.Commons;
	using MGroup.LinearAlgebra.Vectors;

	public interface IGaussSeidelOperable : IIndexable2D
	{
		void GaussSeidelBackwardIteration(IVectorView rhsVector, IVector lhsVector)
		{
			IVector x = lhsVector;
			IVectorView b = rhsVector;
			Preconditions.CheckSquareLinearSystemDimensions(this, x, b);
			int n = this.NumRows;
			for (int i = n - 1; i >= 0; --i)
			{
				double sum = b[i];
				int j;
				for (j = n - 1; j > i; --j)
				{
					sum -= this[i, j] * x[j];
				}

				double diagEntry = this[i, i];

				for (j = i - 1; j >= 0; --j)
				{
					sum -= this[i, j] * x[j];
				}

				x.Set(i, sum / diagEntry);
			}
		}

		void GaussSeidelForwardIteration(IVectorView rhsVector, IVector lhsVector)
		{
			IVector x = lhsVector;
			IVectorView b = rhsVector;
			Preconditions.CheckSquareLinearSystemDimensions(this, x, b);
			int n = this.NumRows;
			for (int i = 0; i < n; ++i)
			{
				double sum = b[i];
				int j;
				for (j = 0; j < i; ++j)
				{
					sum -= this[i, j] * x[j];
				}

				double diagEntry = this[i, i];

				for (j = i + 1; j < n; ++j)
				{
					sum -= this[i, j] * x[j];
				}

				x.Set(i, sum / diagEntry);
			}
		}
	}
}
