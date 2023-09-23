namespace MGroup.LinearAlgebra.Iterative.Termination
{
	using System;
	using System.Collections.Generic;
	using System.Text;

	using MGroup.LinearAlgebra.Vectors;

	public interface ISolutionConvergenceCriterion
	{
		double CalculateConvergenceMetric(IVectorView currentSolution, IVectorView previousSolution);
		string DescribeConvergenceCriterion(double tolerance);
	}
}
