namespace MGroup.LinearAlgebra.Iterative.Preconditioning.Stationary
{
	using System;

	using MGroup.LinearAlgebra.Matrices;

	/// <summary>
	/// Applies the following preconditioning: [1/(w*(2-w) * (D+w*L) * inv(D) * (D+w*U)] * x = y ,
	/// where x is the unknown vector and w the relaxation factor (scalar).
	/// </summary>
	public class SsorPreconditionerCsr : CsrStationaryPreconditionerBase
	{
		private readonly double relaxationFactor;
		private double[] workArray;

		/// <summary>
		/// Creates a new <see cref="SsorPreconditionerCsr"/> with the specified settings.
		/// </summary>
		/// <param name="relaxationFactor">
		/// The w factor used in the SSOR preconditioning: [1/(w*(2-w) * (D+w*L) * inv(D) * (D+w*U)] * x = y.
		/// There is an optimum value for each linear system, but it is unknown apriori. Usually values in [1.0, 2.0) are used.
		/// </param>
		/// <param name="numApplications">
		/// How many times to apply the preconditioning at each iteration of the external iterative algorithm.
		/// </param>
		public SsorPreconditionerCsr(double relaxationFactor, int numApplications = 1)
			: base(numApplications)
		{
			this.relaxationFactor = relaxationFactor;
		}

		public override IPreconditioner CopyWithInitialSettings() => new SsorPreconditionerCsr(relaxationFactor);

		public override void UpdateMatrix(IReadOnlyMatrix matrix, bool isPatternModified)
		{
			base.UpdateMatrix(matrix, isPatternModified);
			if (isPatternModified || workArray == null)
			{
				workArray = new double[matrix.NumColumns];
			}
		}

		protected override void SolveLinearSystemInternal(double[] rhsVector, double[] lhsVector)
		{
			provider.CsrSsorPrecond(matrix.NumRows, matrix.RawValues, matrix.RawRowOffsets, matrix.RawColIndices,
				diagonalOffsets, rhsVector, lhsVector, relaxationFactor, workArray);
			//provider.CsrSorForwardPrecond(matrix.NumRows, matrix.RawValues, matrix.RawRowOffsets, matrix.RawColIndices,
			//	diagonalOffsets, rhsVector, lhsVector, relaxationFactor);
			//provider.CsrSorBackPrecond(matrix.NumRows, matrix.RawValues, matrix.RawRowOffsets, matrix.RawColIndices,
			//	diagonalOffsets, rhsVector, lhsVector, relaxationFactor);
		}
	}
}
