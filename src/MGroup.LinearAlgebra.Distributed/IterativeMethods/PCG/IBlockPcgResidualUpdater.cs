using MGroup.MSolve.Solution.LinearSystem;

namespace MGroup.LinearAlgebra.Distributed.IterativeMethods.PCG
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
        void UpdateResidual(BlockPcgAlgorithm pcg, IGlobalVector residual);
    }
}
