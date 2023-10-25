namespace MGroup.LinearAlgebra.Iterative.AlgebraicMultiGrid
{
	using System;
	using System.Collections.Generic;
	using System.Runtime.CompilerServices;
	using System.Text;

	using DotNumerics.Optimization.LBFGSB;

	using MGroup.LinearAlgebra.Iterative.AlgebraicMultiGrid.Smoothing;
	using MGroup.LinearAlgebra.Iterative.GaussSeidel;
	using MGroup.LinearAlgebra.Matrices;
	using MGroup.LinearAlgebra.Vectors;

	public class GaussSeidelSmoother : IMultigridSmoother
	{
		protected readonly IGaussSeidelIteration _gaussSeidelIteration;
		private readonly GaussSeidelSweepDirection _sweep;
		private readonly int _numIterations;
		private bool _disposed;


		public GaussSeidelSmoother(IGaussSeidelIteration gaussSeidelIteration, GaussSeidelSweepDirection sweep, int numIterations)
		{
			this._gaussSeidelIteration = gaussSeidelIteration;
			this._sweep = sweep;
			this._numIterations = numIterations;
			this._disposed = false;
		}

		public void Dispose()
		{
			if (!_disposed)
			{
				_gaussSeidelIteration.Dispose();
				_disposed = true;
			}
		}

		public void Initialize(IMatrixView matrix)
		{
			CheckDisposed();
			_gaussSeidelIteration.Initialize(matrix);
		}

		public void Smooth(IVectorView rhs, IVector lhs)
		{
			CheckDisposed();
			for (int t = 0; t < _numIterations; t++) 
			{
				_sweep.Apply(_gaussSeidelIteration, rhs, lhs);
			}
		}

		[MethodImpl(MethodImplOptions.AggressiveInlining)]
		private void CheckDisposed()
		{
			if (_disposed)
			{
				throw new ObjectDisposedException(this.GetType().Name);
			}
		}

		public class Builder : IMultigridSmootherBuilder
		{
			private readonly IGaussSeidelIterationBuilder _gaussSeidelIterationBuilder;
			private readonly GaussSeidelSweepDirection _sweep;
			private readonly int _numIterations;

			public Builder(IGaussSeidelIterationBuilder gaussSeidelIterationBuilder, GaussSeidelSweepDirection sweep, 
				int numIterations)
			{
				this._gaussSeidelIterationBuilder = gaussSeidelIterationBuilder;
				this._sweep = sweep;
				this._numIterations = numIterations;
			}

			public IMultigridSmoother Create()
			{
				IGaussSeidelIteration gsIter = _gaussSeidelIterationBuilder.Create();
				return new GaussSeidelSmoother(gsIter, _sweep, _numIterations);
			}
		}
	}
}
