using System;
using System.Collections.Generic;
using System.Text;

namespace MGroup.LinearAlgebra.Distributed.IterativeMethods.PCG
{
    public class PureResidualConvergence : IPcgResidualConvergence
    {
        private double denominator;

        public double EstimateResidualNormRatio(PcgAlgorithmBase pcg) => pcg.Residual.Norm2() / denominator;

        public void Initialize(PcgAlgorithmBase pcg) => denominator = pcg.Rhs.Norm2();
    }
}
