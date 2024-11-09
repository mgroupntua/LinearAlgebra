namespace MGroup.LinearAlgebra.Iterative.Preconditioning.Stationary
{
	using MGroup.LinearAlgebra.Commons;
	using MGroup.LinearAlgebra.Exceptions;
	using MGroup.LinearAlgebra.Matrices;
	using MGroup.LinearAlgebra.Implementations.Managed;
	using MGroup.LinearAlgebra.Vectors;

	public abstract class CsrStationaryPreconditionerBase : IPreconditioner
	{
		protected readonly int numApplications;
		protected readonly StationaryIterationManagedProvider provider;
		protected CsrMatrix matrix;
		protected int[] diagonalOffsets;

		public CsrStationaryPreconditionerBase(int numApplications)
		{
			this.numApplications = numApplications;
			provider = new StationaryIterationManagedProvider();
		}

		public virtual void UpdateMatrix(IMatrixView matrix, bool isPatternModified)
		{
			Preconditions.CheckSquare(matrix);
			if (matrix is CsrMatrix csrMatrix)
			{
				this.matrix = csrMatrix;

				if (isPatternModified || diagonalOffsets == null)
				{
					diagonalOffsets = provider.LocateDiagonalOffsetsCsr(
						matrix.NumRows, csrMatrix.RawRowOffsets, csrMatrix.RawColIndices);
				}
			}
			else
			{
				throw new InvalidSparsityPatternException(GetType().Name + " can be used only for matrices in CSR format.");
			}
		}

		public void SolveLinearSystem(IVectorView rhsVector, IVector lhsVector)
		{
			if (rhsVector is Vector rhs && lhsVector is Vector lhs)
			{
				for (var i = 0; i < numApplications; ++i)
				{
					SolveLinearSystemInternal(rhs.RawData, lhs.RawData);
				}
			}
			else
			{
				throw new NonMatchingFormatException("Input and output vectors must be dense.");
			}
		}

		public abstract IPreconditioner CopyWithInitialSettings();

		protected abstract void SolveLinearSystemInternal(double[] rhsVector, double[] lhsVector);
	}
}
