namespace MGroup.LinearAlgebra.Iterative.PreconditionedConjugateGradient.Logging
{
	using System.Text;

	/// <summary>
	/// Manages multiple <see cref="IPcgLogger"/> applied simultaneously.
	/// </summary>
	public class CompositePcgLogger : IPcgLogger
	{
		private readonly IPcgLogger[] loggers;

		public CompositePcgLogger(params IPcgLogger[] loggers)
		{
			this.loggers = loggers;
		}

		public void Clear()
		{
			foreach (IPcgLogger logger in loggers)
			{
				logger.Clear();
			}
		}

		public void Log(PcgAlgorithmBase pcg)
		{
			foreach (IPcgLogger logger in loggers)
			{
				logger.Log(pcg);
			}
		}

		public string Report()
		{
			var output = new StringBuilder();
			foreach (IPcgLogger logger in loggers)
			{
				output.AppendLine();
				output.AppendLine($"*** Report from {logger.GetType()} ***");
				output.AppendLine(logger.Report());
			}

			return output.ToString();
		}
	}
}
