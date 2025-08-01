namespace MGroup.LinearAlgebra.Iterative.PreconditionedConjugateGradient
{
	using MGroup.LinearAlgebra.Vectors;

	/// <summary>
	/// Updates the residual vector r using either the standard operation r = r - α * A*d or the 
	/// exact formula r = b - A*x.
	/// </summary>
	public interface IPcgResidualUpdater : ISettingsCopiable<IPcgResidualUpdater>
	{
		/// <summary>
		/// Update the residual vector r.
		/// </summary>
		/// <param name="pcg">The Preconditioned Conjugate Gradient algorithm that uses this object.</param>
		/// <param name="residual">The current residual vector r to modify.</param>
		void UpdateResidual(PcgAlgorithmBase pcg, IVector residual);
	}
}
