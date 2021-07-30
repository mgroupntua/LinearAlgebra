using System;
using System.Collections.Generic;
using System.Text;

namespace MGroup.LinearAlgebra.Iterative.Termination
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
            this.residualDotProductsHistory = new List<double>();
        }

        public bool HasStagnated()
        {
            int numIterations = residualDotProductsHistory.Count;
            if (numIterations < iterationSpan) return false; // Not enough data yet
            double oldError = residualDotProductsHistory[numIterations - iterationSpan];
            double newError = residualDotProductsHistory[numIterations - 1];
            double relativeReduction = (oldError - newError) / oldError;
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
            int t = 0;
            while (t < iterationSpan)
            {
                double current = residualDotProductsHistory[t];
                double next = residualDotProductsHistory[t + 1];
                double reduction = (current - next) / current;
                if (reduction > 0) return reduction;
                else ++t;
            }
            // At this point PCG has made no improvement. Thus it diverges.
            throw new Exception("PCG diverges");
        }
    }
}
