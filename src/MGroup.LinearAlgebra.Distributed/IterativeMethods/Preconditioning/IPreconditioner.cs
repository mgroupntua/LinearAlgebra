using System;
using System.Collections.Generic;
using System.Text;

using MGroup.MSolve.Solution.LinearSystem;

//TODOMPI: Make it generic, fix comments. IIndexable1D is no longer valid for iterative method vectors, since random access does 
//      not make sense for distributed vectors. IBounded1D would be better here, but I do not know if forcing distributed vectors
//      to know the whole length is a good design. Also fix XML comments.
namespace MGroup.LinearAlgebra.Distributed.IterativeMethods.Preconditioning
{
    public interface IPreconditioner
    {
        /// <summary>
        /// Applies the preconditioner by solving the system: M * v = w during the preconditioning step of an iterative 
        /// algorithm, where M is the preconditioner and the definition of the vectors v, w depends on the iterative algorithm. 
        /// This is equivalent to evaluating: inverse(M) * w.
        /// </summary>
        /// <param name="input">The right hand side vector of the system M * v = w. It will not be modified.</param>
        /// <param name="output">
        /// The left hand side vector of the system M * v = w. It will be overwritten by the solution of the linear system.
        /// </param>
        /// <exception cref="LinearAlgebraExtensions.NonMatchingDimensionsException">
        /// Thrown if the <see cref="IIndexable1D.Length"/> of <paramref name="input"/> or <paramref name="output"/> 
        /// is different than the number of rows of this <see cref="IPreconditioner"/>.
        /// </exception>
        void Apply(IGlobalVector input, IGlobalVector output);
    }
}
