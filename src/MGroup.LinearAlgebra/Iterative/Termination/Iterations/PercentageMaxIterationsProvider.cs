﻿using System;

using MGroup.LinearAlgebra.Commons;

namespace MGroup.LinearAlgebra.Iterative.Termination.Iterations
{
	/// <summary>
	/// A <see cref="IMaxIterationsProvider"/> implementation that will use a percentage of the matrix order as the max number
	/// of iterations.
	/// Authors: Serafeim Bakalakos
	/// </summary>
	public class PercentageMaxIterationsProvider : IMaxIterationsProvider
	{
		private readonly double maxIterationsOverMatrixOrder;

		/// <summary>
		/// Initializes a new instance of <see cref="PercentageMaxIterationsProvider"/> with the specified settings.
		/// </summary>
		/// <param name="maxIterationsOverMatrixOrder">
		/// The percentage of the matrix order that will be set as max iterations. It will be rounded up. Constraints:
		/// 0.0 &lt; <paramref name="maxIterationsOverMatrixOrder"/> &lt;= 1.0.
		/// </param>
		public PercentageMaxIterationsProvider(double maxIterationsOverMatrixOrder)
		{
			if (maxIterationsOverMatrixOrder <= 0.0 || maxIterationsOverMatrixOrder > 1.0) throw new ArgumentException(
				"The ratio of max iterations / matrix order must belong to the interval (0.0, 1.0],"
				+ $" but was {maxIterationsOverMatrixOrder}");
			this.maxIterationsOverMatrixOrder = maxIterationsOverMatrixOrder;
		}

		/// <summary>
		/// See <see cref="IMaxIterationsProvider.GetMaxIterations(int)"/>
		/// </summary>
		public int GetMaxIterations(int matrixOrder)
		{
			if (maxIterationsOverMatrixOrder == 1.0) return matrixOrder;
			else return (int)Math.Ceiling(maxIterationsOverMatrixOrder * matrixOrder);
		}
	}
}
