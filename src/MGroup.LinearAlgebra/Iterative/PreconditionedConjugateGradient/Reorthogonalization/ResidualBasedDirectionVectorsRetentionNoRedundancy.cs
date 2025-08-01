namespace MGroup.LinearAlgebra.Iterative.PreconditionedConjugateGradient.Reorthogonalization
{
	public class ResidualBasedDirectionVectorsRetentionNoRedundancy : IDirectionVectorsRetention
	{
		private readonly double minResidualNormRatioToKeep;

		private ReorthogonalizedPcg pcg;

		public ResidualBasedDirectionVectorsRetentionNoRedundancy(double minResidualNormRatioToKeep)
		{
			this.minResidualNormRatioToKeep = minResidualNormRatioToKeep;
		}

		public IDirectionVectorsRetention CopyWithInitialSettings()
			=> new ResidualBasedDirectionVectorsRetentionNoRedundancy(minResidualNormRatioToKeep);

		public void DiscardDirectionVectors() { }

		public void Initialize(ReorthogonalizedPcg pcg)
		{
			this.pcg = pcg;
		}

		public bool KeepUsingReorthogonalization()
		{
			if (double.IsNaN(pcg.ResidualNormRatio))
			{
				return true; // First iteration - not initialized yet.
			}
			else
			{
				return pcg.ResidualNormRatio > minResidualNormRatioToKeep;
			}
		}
	}
}
