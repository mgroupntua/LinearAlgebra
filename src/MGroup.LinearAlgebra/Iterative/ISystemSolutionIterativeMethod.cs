namespace MGroup.LinearAlgebra.Iterative
{
	using System;
	using System.Collections.Generic;
	using System.Text;

	using MGroup.LinearAlgebra.Iterative.Preconditioning;
	using MGroup.LinearAlgebra.Matrices;
	using MGroup.LinearAlgebra.Vectors;

	public interface ISystemSolutionIterativeMethod
	{
		void Clear();

		IterativeStatistics Solve(ILinearTransformation matrix, IPreconditioner preconditioner, IVectorView rhs,
			IVector solution, bool initialGuessIsZero);
	}
}
