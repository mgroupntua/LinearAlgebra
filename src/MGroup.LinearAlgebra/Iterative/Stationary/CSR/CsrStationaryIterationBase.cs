namespace MGroup.LinearAlgebra.Iterative.Stationary.CSR
{
	using System.Collections.Generic;

	using MGroup.LinearAlgebra.Commons;
	using MGroup.LinearAlgebra.Exceptions;
	using MGroup.LinearAlgebra.Matrices;
	using MGroup.LinearAlgebra.Implementations.Managed;
	using MGroup.LinearAlgebra.Vectors;

	public abstract class CsrStationaryIterationBase : IStationaryIteration
	{
		protected readonly List<IStationaryIteration> linkedIterations = new List<IStationaryIteration>();
		protected readonly StationaryIterationManagedProvider provider;
		protected CsrMatrix matrix;
		protected int[] diagonalOffsets;

		public CsrStationaryIterationBase()
		{
			provider = new StationaryIterationManagedProvider();
		}

		public void LinkWith(IStationaryIteration other)
		{
			linkedIterations.Add(other);
		}

		public virtual void UpdateMatrix(IMatrixView matrix, bool isPatternModified)
		{
			Preconditions.CheckSquare(matrix);
			if (matrix is CsrMatrix csrMatrix)
			{
				this.matrix = csrMatrix;

				if (isPatternModified || diagonalOffsets == null)
				{
					foreach (IStationaryIteration iteration in linkedIterations)
					{
						if (iteration is CsrStationaryIterationBase csrIteration)
						{
							diagonalOffsets = csrIteration.diagonalOffsets;
							return;
						}
					}

					if (diagonalOffsets == null)
					{
						diagonalOffsets = provider.LocateDiagonalOffsetsCsr(
							matrix.NumRows, csrMatrix.RawRowOffsets, csrMatrix.RawColIndices);
					}
				}
			}
			else
			{
				throw new InvalidSparsityPatternException(GetType().Name + " can be used only for matrices in CSR format.");
			}
		}

		public abstract string Name { get; }

		public abstract IStationaryIteration CopyWithInitialSettings();

		public abstract void Execute(Vector rhs, Vector solution);
	}
}
