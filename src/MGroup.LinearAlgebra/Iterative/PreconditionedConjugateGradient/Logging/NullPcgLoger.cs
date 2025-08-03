namespace MGroup.LinearAlgebra.Iterative.PreconditionedConjugateGradient.Logging
{
	/// <summary>
	/// Does not log anything and does minimal work. To be used as the default <see cref="IPcgLogger"/>,
	/// when another implementation is not common.
	/// </summary>
	public class NullPcgLoger : IPcgLogger
	{
		public void Clear() { }

		public void Log(PcgAlgorithmBase pcg) { }

		public string Report() => string.Empty;
	}
}
