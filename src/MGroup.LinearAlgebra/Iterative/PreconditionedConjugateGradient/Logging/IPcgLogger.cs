
namespace MGroup.LinearAlgebra.Iterative.PreconditionedConjugateGradient.Logging
{
	public interface IPcgLogger
	{
		void Clear();

		void Log(PcgAlgorithmBase pcg);

		string Report();
	}
}
