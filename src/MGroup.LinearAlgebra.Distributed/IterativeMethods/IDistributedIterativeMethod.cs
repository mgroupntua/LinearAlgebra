using System;
using System.Collections.Generic;
using System.Text;
using MGroup.LinearAlgebra.Iterative;
using MGroup.LinearAlgebra.Distributed.IterativeMethods.Preconditioning;
using MGroup.LinearAlgebra.Vectors;

namespace MGroup.LinearAlgebra.Distributed.IterativeMethods
{
	public interface IDistributedIterativeMethod //TODO: Use the regular iterative methods, instead of these classes
	{
		void Clear();

		IterativeStatistics Solve(ILinearTransformation matrix, IPreconditioner preconditioner,
			IVector rhs, IVector solution, bool initialGuessIsZero);
	}
}
