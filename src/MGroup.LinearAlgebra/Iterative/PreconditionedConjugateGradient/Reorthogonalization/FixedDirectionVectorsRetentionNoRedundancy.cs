namespace MGroup.LinearAlgebra.Iterative.PreconditionedConjugateGradient.Reorthogonalization
{
	public class FixedDirectionVectorsRetentionNoRedundancy : IDirectionVectorsRetention
	{
		private readonly int numVectorsToKeep;

		private ReorthogonalizedPcg pcg;

		public FixedDirectionVectorsRetentionNoRedundancy(int numVectorsToKeep)
		{
			this.numVectorsToKeep = numVectorsToKeep;
		}

		public IDirectionVectorsRetention CopyWithInitialSettings()
			=> new FixedDirectionVectorsRetentionNoRedundancy(numVectorsToKeep);

		public void DiscardDirectionVectors() { }

		public void Initialize(ReorthogonalizedPcg pcg)
		{
			this.pcg = pcg;
		}

		public bool KeepUsingReorthogonalization() => pcg.ReorthoCache.Directions.Count < numVectorsToKeep;
	}
}
