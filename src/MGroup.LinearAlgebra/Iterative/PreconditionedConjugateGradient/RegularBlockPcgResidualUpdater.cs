using System;
using System.Collections.Generic;
using System.Text;
using MGroup.LinearAlgebra.Vectors;

//TODO: Duplication between this, the CG and the PCG version
namespace MGroup.LinearAlgebra.Iterative.PreconditionedConjugateGradient
{
    /// <summary>
    /// Updates the residual vector according to the usual CG formula r = r - α * A*d. No corrections are applied.
    /// </summary>
    public class RegularBlockPcgResidualUpdater : IBlockPcgResidualUpdater
    {
        /// <summary>
        /// See <see cref="IBlockPcgResidualUpdater.UpdateResidual(BlockPcgAlgorithm, IVector)"/>
        /// </summary>
        public void UpdateResidual(BlockPcgAlgorithm pcg, IVector residual)
        {
            // Normally the residual vector is updated as: r = r - α * A*d
			residual.CopyFrom(pcg.ResidualOperator.EvaluateVector(pcg.ResidualKernels, pcg.DirectionKernels));  // It didn't multiplied with M, because it shouldn't be
		}

		/// <summary>
		/// See <see cref="IPcgResidualUpdater.UpdateResidual(PcgAlgorithmBase, IVector)"/>
		/// </summary>
		public void UpdateResidual(PcgAlgorithmBase pcg, IVector residual)
		{
			// Normally the residual vector is updated as: r = r - α * A*d
			residual.AxpyIntoThis(pcg.MatrixTimesDirection, -pcg.StepSize);
		}
	}
}
