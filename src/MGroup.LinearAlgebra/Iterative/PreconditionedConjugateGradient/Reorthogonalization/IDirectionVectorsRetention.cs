namespace MGroup.LinearAlgebra.Iterative.PreconditionedConjugateGradient.Reorthogonalization
{
	public interface IDirectionVectorsRetention : ISettingsCopiable<IDirectionVectorsRetention>
	{
		void DiscardDirectionVectors();

		void Initialize(ReorthogonalizedPcg pcg);

		bool KeepUsingReorthogonalization();
	}
}
