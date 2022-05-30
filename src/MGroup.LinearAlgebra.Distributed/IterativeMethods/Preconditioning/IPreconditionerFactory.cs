using System;
using System.Collections.Generic;
using System.Text;

using MGroup.MSolve.Solution.LinearSystem;

//TODO: working with the matrix itself is not always convenient, especially when abstracted behind IMatrixView. E.g. Jacobi
//      preconditioner needs to access the diagonal which is more efficient with the DOK. IncompleteCholesky may also be more 
//      effient on other matrix storage formats than the CSR that will be used for multiplications.
//TODO: This and its implementation may need to be removed, as the creation of preconditioners is extremely problem specific and 
//      there is no benefit in having a general interface, which implies that the preconditioner will be created as a function 
//      of the matrix itself.
namespace MGroup.LinearAlgebra.Distributed.IterativeMethods.Preconditioning
{
    /// <summary>
    /// Classes implementing this interface are responsible for the creation of <see cref="IPreconditioner"/> instances.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public interface IPreconditionerFactory
    {
        /// <summary>
        /// Initializes and returns an <see cref="IPreconditioner"/> for the provided <paramref name="matrix"/>.
        /// </summary>
        /// <param name="matrix">The original matrix, whose preconditioner will be created.</param>
        IPreconditioner CreatePreconditionerFor(ILinearTransformation matrix);
    }
}
