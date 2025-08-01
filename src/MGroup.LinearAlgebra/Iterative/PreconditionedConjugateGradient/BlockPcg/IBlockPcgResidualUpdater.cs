namespace MGroup.LinearAlgebra.Iterative.PreconditionedConjugateGradient.BlockPcg
{
	using MGroup.LinearAlgebra.Iterative;
	using MGroup.LinearAlgebra.Vectors;

	/// <summary>
	/// Abstraction to update the residual vector r.
	/// </summary>
	public interface IBlockPcgResidualUpdater : ISettingsCopiable<IBlockPcgResidualUpdater>
	{
		/// <summary>
		/// Update the residual vector r.
		/// </summary>
		/// <param name="pcg">The Block Preconditioned Conjugate Gradient algorithm that uses this object.</param>
		/// <param name="residual">The current residual vector r to modify.</param>
		void UpdateResidual(BlockPcgAlgorithm pcg, IVector residual);
	}
}
