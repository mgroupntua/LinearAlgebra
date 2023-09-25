namespace MGroup.LinearAlgebra.Iterative.AlgebraicMultiGrid
{
	using System;
	using System.Collections.Generic;
	using System.Text;

	using MGroup.LinearAlgebra.Iterative.AlgebraicMultiGrid.Smoothing;
	using MGroup.LinearAlgebra.Iterative.GaussSeidel;

	public class PodAmg
	{
		private Func<IMultigridSmoother> choosePresmoother;
		private Func<IMultigridSmoother> choosePostsmoother;

	}
}
