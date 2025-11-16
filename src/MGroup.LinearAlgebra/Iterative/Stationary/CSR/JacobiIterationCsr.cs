namespace MGroup.LinearAlgebra.Iterative.Stationary.CSR
{
	using System;
	using System.Collections.Generic;
	using System.Text;

	using MGroup.LinearAlgebra.Commons;
	using MGroup.LinearAlgebra.Exceptions;
	using MGroup.LinearAlgebra.Matrices;
	using MGroup.LinearAlgebra.Vectors;

	public class JacobiIterationCsr : CsrStationaryIterationBase
	{
		private double[] workArray;

		public override string Name => "Jacobi";

		public override IStationaryIteration CopyWithInitialSettings() => new JacobiIterationCsr();

		public override void Execute(Vector rhs, Vector solution)
		{
			provider.CsrJacobi(matrix.NumRows, matrix.RawValues, matrix.RawRowOffsets, matrix.RawColIndices,
					diagonalOffsets, rhs.RawData, solution.RawData, workArray);
		}

		public override void UpdateMatrix(IReadOnlyMatrix matrix, bool isPatternModified)
		{
			base.UpdateMatrix(matrix, isPatternModified);
			if (isPatternModified || workArray == null)
			{
				workArray = new double[matrix.NumColumns];
			}
		}
	}
}
