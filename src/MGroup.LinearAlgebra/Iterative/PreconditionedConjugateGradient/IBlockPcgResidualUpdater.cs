using System;
using System.Collections.Generic;
using System.Text;
using MGroup.LinearAlgebra.Vectors;

namespace MGroup.LinearAlgebra.Iterative.PreconditionedConjugateGradient
{
    /// <summary>
    /// Abstraction to update the residual vector r 
    /// </summary>
    public interface IBlockPcgResidualUpdater : IPcgResidualUpdater
    {
        /// <summary>
        /// Update the residual vector r.
        /// </summary>
        ///<param name="pcg">The Block Preconditioned Conjugate Gradient algorithm that uses this object.</param>
        /// <param name="residual">The current residual vector r to modify.</param>
        void UpdateResidual(BlockPcgAlgorithm pcg, IVector residual);
    }
}
