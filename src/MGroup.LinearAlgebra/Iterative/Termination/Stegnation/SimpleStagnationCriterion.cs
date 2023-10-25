using System;
using System.Collections.Generic;
using System.Text;

namespace MGroup.LinearAlgebra.Iterative.Termination.Stegnation
{
	public class SimpleStagnationCriterion : IStagnationCriterion
	{
		private readonly int iterationSpan;
		private double relativeImprovementTolerance;
		private List<double> residualDotProductsHistory;

		public SimpleStagnationCriterion(int iterationSpan, double relativeImprovementTolerance = -1)
		{
			this.iterationSpan = iterationSpan;
			this.relativeImprovementTolerance = relativeImprovementTolerance;
			residualDotProductsHistory = new List<double>();
		}

		public bool HasStagnated()
		{
			var numIterations = residualDotProductsHistory.Count;
			if (numIterations < iterationSpan) return false; // Not enough data yet
			var oldError = residualDotProductsHistory[numIterations - iterationSpan];
			var newError = residualDotProductsHistory[numIterations - 1];
			var relativeReduction = (oldError - newError) / oldError;
			if (relativeImprovementTolerance == -1)
			{
				relativeImprovementTolerance = 1E-3 * CalcInitialErrorReduction();
			}
			if (relativeReduction <= relativeImprovementTolerance) return true;
			else return false;
		}

		public void StoreInitialError(double initialError)
		{
			residualDotProductsHistory.Clear();
			residualDotProductsHistory.Add(initialError);
		}

		public void StoreNewError(double currentError)
		{
			residualDotProductsHistory.Add(currentError);
		}

		private double CalcInitialErrorReduction()
		{
			var t = 0;
			while (t < iterationSpan)
			{
				var current = residualDotProductsHistory[t];
				var next = residualDotProductsHistory[t + 1];
				var reduction = (current - next) / current;
				if (reduction > 0) return reduction;
				else ++t;
			}
			// At this point PCG has made no improvement. Thus it diverges.
			throw new Exception("PCG diverges");
		}
	}
}
