namespace MGroup.LinearAlgebra.Iterative.PreconditionedConjugateGradient
{
	public class PureResidualConvergence : IPcgResidualConvergence
	{
		private double denominator;

		public IPcgResidualConvergence CopyWithInitialSettings()
		{
			var clone = new PureResidualConvergence();
			clone.denominator = this.denominator;
			return clone;
		}

		public double EstimateResidualNormRatio(PcgAlgorithmBase pcg) => pcg.Residual.Norm2() / denominator;

		public void Initialize(PcgAlgorithmBase pcg) => denominator = pcg.Rhs.Norm2();
	}
}
