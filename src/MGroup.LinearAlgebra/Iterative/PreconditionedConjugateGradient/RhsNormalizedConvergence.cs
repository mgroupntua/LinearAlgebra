namespace MGroup.LinearAlgebra.Iterative.PreconditionedConjugateGradient
{
	using System;

	public class RhsNormalizedConvergence : IPcgResidualConvergence
	{
		private double denominator;

		public IPcgResidualConvergence CopyWithInitialSettings()
		{
			var clone = new RhsNormalizedConvergence();
			clone.denominator = this.denominator;
			return clone;
		}

		public double EstimateResidualNormRatio(PcgAlgorithmBase pcg) => Math.Sqrt(pcg.ResDotPrecondRes) / denominator;

		public void Initialize(PcgAlgorithmBase pcg) => denominator = Math.Sqrt(pcg.Rhs.Norm2());
	}
}
