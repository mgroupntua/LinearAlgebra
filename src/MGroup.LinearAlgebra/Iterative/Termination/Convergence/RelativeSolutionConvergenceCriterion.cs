namespace MGroup.LinearAlgebra.Iterative.Termination.Convergence
{
	using System;
	using System.Collections.Generic;
	using System.Text;

	using MGroup.LinearAlgebra.Vectors;

	/// <summary>
	/// Convergence criterion: norm2(x - x_previous) / norm2(x) &lt; tolerance, where x is the solution vector.
	/// </summary>
	public class RelativeSolutionConvergenceCriterion : ISolutionConvergenceCriterion
	{
		public double CalculateConvergenceMetric(IVectorView currentSolution, IVectorView previousSolution)
		{
			//TODO: The next can be optimized to not create a new vector (using SubtractIntoThis) in some cases.
			//		E.g. in Gauss-Seidel the previousSolution vector is no longer necessary and can be overwritten.
			var num = previousSolution.Subtract(currentSolution).Norm2();
			var den = currentSolution.Norm2();
			if (den != 0)
			{
				return num / den;
			}
			else
			{
				return num;
			}
		}

		public string DescribeConvergenceCriterion(double convergenceTolerance) => "norm2(x - x_previous) / norm2(x) < " + convergenceTolerance;
	}
}
