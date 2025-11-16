namespace MGroup.LinearAlgebra.Iterative.Termination.Convergence
{
	using MGroup.LinearAlgebra.Vectors;

	/// <summary>
	/// Useful when Gauss-Seidel needs to run for a specified number of iterations, without converging. 
	/// This class avoids redundant operations to calculate the convergence metric at each iteration.
	/// </summary>
	public class SolutionNeverConvergesCriterion : ISolutionConvergenceCriterion
	{
		public double CalculateConvergenceMetric(IReadOnlyVector currentSolution, IReadOnlyVector previousSolution) => double.MaxValue;

		public ISolutionConvergenceCriterion CopyWithInitialSettings() => new SolutionNeverConvergesCriterion();

		public string DescribeConvergenceCriterion(double tolerance)
			=> "No convergence criterion specified. Iterative solution algorithm stops, when max iterations are reached.";
	}
}
