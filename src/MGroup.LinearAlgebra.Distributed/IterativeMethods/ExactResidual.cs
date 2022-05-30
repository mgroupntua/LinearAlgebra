using MGroup.LinearAlgebra.Vectors;
using MGroup.MSolve.Solution.LinearSystem;

namespace MGroup.LinearAlgebra.Distributed.IterativeMethods
{
    internal static class ExactResidual
    {
        internal static IGlobalVector Calculate(ILinearTransformation matrix,
            IGlobalVector rhs, IGlobalVector solution)
        {
            IGlobalVector residual = rhs.CreateZero();
            Calculate(matrix, rhs, solution, residual);
            return residual;
        }

        internal static void Calculate(ILinearTransformation matrix, IGlobalVector rhs,
            IGlobalVector solution, IGlobalVector residual)
        {
            //TODO: There is a BLAS operation y = y + a * A*x, that would be perfect for here. rhs.Copy() and then that.
            matrix.MultiplyVector(solution, residual);
            residual.LinearCombinationIntoThis(-1.0, rhs, 1.0);
        }

    }
}
