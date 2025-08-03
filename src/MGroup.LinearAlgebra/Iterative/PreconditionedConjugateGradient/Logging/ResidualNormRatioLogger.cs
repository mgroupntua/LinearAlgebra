namespace MGroup.LinearAlgebra.Iterative.PreconditionedConjugateGradient.Logging
{
	using System.Collections.Generic;
	using System.Text;

	/// <summary>
	/// Logs the normalized value of ||r^T * inv(M) * r||, where r = residual and M = preconditioner matrix. 
	/// </summary>
	public class ResidualNormRatioLogger : IPcgLogger
	{
		private readonly List<double> residualNormRatios = [];

		public void Clear()
		{
			residualNormRatios.Clear();
		}

		public void Log(PcgAlgorithmBase pcg)
		{
			residualNormRatios.Add(pcg.ConvergenceStrategy.EstimateResidualNormRatio(pcg));
		}

		public string Report()
		{
			var output = new StringBuilder();
			output.AppendLine("Normalized values of r^T * inv(M) * r, where r = residual and M = preconditioner matrix:");
			for (int t = 0; t < residualNormRatios.Count; t++)
			{
				output.AppendLine($"Iteration {t}: {residualNormRatios[t]}");
			}

			return output.ToString();
		}
	}
}
