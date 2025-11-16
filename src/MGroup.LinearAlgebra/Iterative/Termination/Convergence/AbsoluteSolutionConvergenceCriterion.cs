namespace MGroup.LinearAlgebra.Iterative.Termination.Convergence
{
	using MGroup.LinearAlgebra.Vectors;

	/// <summary>
	/// Convergence criterion: norm2(x - x_previous) &lt; tolerance, where x is the solution vector.
	/// </summary>
	public class AbsoluteSolutionConvergenceCriterion : ISolutionConvergenceCriterion
	{
		public double CalculateConvergenceMetric(IReadOnlyVector currentSolution, IReadOnlyVector previousSolution)
		{
			//TODO: The next can be optimized to not create a new vector (using SubtractIntoThis) in some cases.
			//		E.g. in Gauss-Seidel the previousSolution vector is no longer necessary and can be overwritten.
			return previousSolution.Subtract(currentSolution).Norm2();
		}

		public ISolutionConvergenceCriterion CopyWithInitialSettings() => new AbsoluteSolutionConvergenceCriterion();

		public string DescribeConvergenceCriterion(double convergenceTolerance) => "norm2(x - x_previous) < " + convergenceTolerance;
	}
}
