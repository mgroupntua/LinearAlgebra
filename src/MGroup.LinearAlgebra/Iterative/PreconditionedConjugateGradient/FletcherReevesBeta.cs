namespace MGroup.LinearAlgebra.Iterative.PreconditionedConjugateGradient
{
	using MGroup.LinearAlgebra.Vectors;

	/// <summary>
	/// Fletcher-Reeves: beta = (rNew * inv(M) * rNew) / (rOld * inv(M) * rOld).
	/// This is the simplest formula to calculate PCG's beta parameter and does not require any extra memory or calculations.
	/// </summary>
	public class FletcherReevesBeta : IPcgBetaParameterCalculation
	{
		/// <summary>
		/// See <see cref="IPcgBetaParameterCalculation.CalculateBeta(PcgAlgorithmBase)"/>.
		/// </summary>
		public double CalculateBeta(PcgAlgorithmBase pcg) => pcg.ResDotPrecondRes / pcg.ResDotPrecondResOld;

		public IPcgBetaParameterCalculation CopyWithInitialSettings() => new FletcherReevesBeta();

		/// <summary>
		/// See <see cref="IPcgBetaParameterCalculation.Initialize(PcgAlgorithmBase)"/>.
		/// </summary>
		public void Initialize(PcgAlgorithmBase pcg) { } // Do nothing
	}
}
