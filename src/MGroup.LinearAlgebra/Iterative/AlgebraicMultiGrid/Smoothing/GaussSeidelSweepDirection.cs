namespace MGroup.LinearAlgebra.Iterative.AlgebraicMultiGrid.Smoothing
{
	using System;
	using System.Collections.Generic;
	using System.Text;

	using MGroup.LinearAlgebra.Iterative.GaussSeidel;
	using MGroup.LinearAlgebra.Vectors;

	public abstract class GaussSeidelSweepDirection
	{
		public static readonly GaussSeidelSweepDirection Backward = new GaussSeidelSweepDirection.BackwardSweep();
		public static readonly GaussSeidelSweepDirection Forward = new GaussSeidelSweepDirection.ForwardSweep();
		public static readonly GaussSeidelSweepDirection Symmetric = new GaussSeidelSweepDirection.SymmetricSweep();

		protected GaussSeidelSweepDirection()
		{
		}

		public abstract void Apply(IGaussSeidelIteration _gsIteration, IVectorView rhs, IVector lhs);


		private class BackwardSweep : GaussSeidelSweepDirection
		{
			public override void Apply(IGaussSeidelIteration _gsIteration, IVectorView rhs, IVector lhs)
			{
				_gsIteration.GaussSeidelBackwardIteration(rhs, lhs);
			}
		}

		private class ForwardSweep : GaussSeidelSweepDirection
		{
			public override void Apply(IGaussSeidelIteration _gsIteration, IVectorView rhs, IVector lhs)
			{
				_gsIteration.GaussSeidelForwardIteration(rhs, lhs);
			}
		}

		private class SymmetricSweep : GaussSeidelSweepDirection
		{
			public override void Apply(IGaussSeidelIteration _gsIteration, IVectorView rhs, IVector lhs)
			{
				_gsIteration.GaussSeidelForwardIteration(rhs, lhs);
				_gsIteration.GaussSeidelBackwardIteration(rhs, lhs);
			}
		}
	}
}
