using System;
using System.Collections.Generic;
using System.Text;

namespace MGroup.LinearAlgebra.Distributed.IterativeMethods.PCG
{
    public class RhsNormalizedConvergence : IPcgResidualConvergence
    {
        private double denominator;

        public double EstimateResidualNormRatio(PcgAlgorithmBase pcg) => Math.Sqrt(pcg.ResDotPrecondRes) / denominator;

        public void Initialize(PcgAlgorithmBase pcg) => denominator = Math.Sqrt(pcg.Rhs.Norm2());
    }
}
