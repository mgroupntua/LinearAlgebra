using System;
using System.Collections.Generic;
using System.Text;

//TODO: Combine this with convergence criterion
namespace MGroup.LinearAlgebra.Iterative.Termination
{
    public interface IStagnationCriterion
    {
        bool HasStagnated();

        void StoreInitialError(double initialError);

        void StoreNewError(double currentError);
    }
}
