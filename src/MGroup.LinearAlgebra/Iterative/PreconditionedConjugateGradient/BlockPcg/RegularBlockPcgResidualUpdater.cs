//TODO: Duplication between this, the CG and the PCG version
namespace MGroup.LinearAlgebra.Iterative.PreconditionedConjugateGradient.BlockPcg
{
	using MGroup.LinearAlgebra.Vectors;

	/// <summary>
	/// Updates the residual vector according to the usual CG formula r = r - α * A*d. No corrections are applied.
	/// </summary>
	public class RegularBlockPcgResidualUpdater : IBlockPcgResidualUpdater
	{
		public IBlockPcgResidualUpdater CopyWithInitialSettings() => new RegularBlockPcgResidualUpdater();

		/// <summary>
		/// See <see cref="IBlockPcgResidualUpdater.UpdateResidual(BlockPcgAlgorithm, IVector)"/>
		/// </summary>
		public void UpdateResidual(BlockPcgAlgorithm pcg, IVector residual)
		{
			// Normally the residual vector is updated as: r = r - α * A*d
			residual.CopyFrom(pcg.ResidualOperator.EvaluateVector(pcg.ResidualKernels, pcg.DirectionKernels));  // It didn't multiplied with M, because it shouldn't be
		}
	}
}
