namespace MGroup.LinearAlgebra.Iterative.PreconditionedConjugateGradient.Reorthogonalization
{
	public class ResidualBasedDirectionVectorsRetention : IDirectionVectorsRetention
	{
		private readonly double minResidualNormRatioToKeep;

		private ReorthogonalizedPcg pcg;
		private int numVectorsToKeep = 0;

		public ResidualBasedDirectionVectorsRetention(double minResidualNormRatioToKeep)
		{
			this.minResidualNormRatioToKeep = minResidualNormRatioToKeep;
		}

		public IDirectionVectorsRetention CopyWithInitialSettings()
			=> new ResidualBasedDirectionVectorsRetention(minResidualNormRatioToKeep);

		public void DiscardDirectionVectors()
		{
			int numVectorsToDiscard = pcg.ReorthoCache.Directions.Count - numVectorsToKeep;
			pcg.ReorthoCache.RemoveNewDirectionVectorData(numVectorsToDiscard);
		}

		public void Initialize(ReorthogonalizedPcg pcg)
		{
			this.pcg = pcg;
		}

		public bool KeepUsingReorthogonalization()
		{
			if (pcg.ResidualNormRatio >= minResidualNormRatioToKeep)
			{
				numVectorsToKeep = pcg.ReorthoCache.Directions.Count;
			}

			return true;
		}
	}
}
