using System;
using System.Collections.Generic;
using System.Text;

namespace MGroup.LinearAlgebra.Iterative.Termination
{
    public class NullStagnationCriterion : IStagnationCriterion
    {
        public bool HasStagnated() => false;

        public void StoreInitialError(double initialError) { }

        public void StoreNewError(double currentError) { }

    }
}
