namespace MGroup.LinearAlgebra.Iterative.Termination
{
	using System;
	using System.Collections.Generic;
	using System.Text;

	using MGroup.LinearAlgebra.Vectors;

	/// <summary>
	/// Useful when Gauss-Seidel needs to run for a specified number of iterations, without converging. 
	/// This class avoids redundant operations to calculate the convergence metric at each iteration.
	/// </summary>
	public class SolutionNeverConvergesCriterion : ISolutionConvergenceCriterion
	{
		public double CalculateConvergenceMetric(IVectorView currentSolution, IVectorView previousSolution) => double.MaxValue;

		public string DescribeConvergenceCriterion(double tolerance) 
			=> "No convergence criterion specified. Iterative solution algorithm stops, when max iterations are reached.";
	}
}
