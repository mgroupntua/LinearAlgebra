namespace MGroup.LinearAlgebra.Iterative
{
	using MGroup.LinearAlgebra.Vectors;

	/** Interface to check if an iterative algorithm converged. */
	public interface IResidualConvergence
	{
		
		/** Checks if iterative algorithm converged.
			* <param name="rr">The residual^2</param>
			* <returns>True if converged, false otherwise</returns>
			*/
		public bool isConverged(double rr);

		/** Checks if iterative algorithm converged.
			* Better avoid to use this method, but it exists if you want fancy residual checks.
			* <param name="r">The residual vector</param>
			* <returns>True if converged, false otherwise</returns>
			*/
		public bool isConverged(IVectorView r) { return isConverged(r.DotProduct(r)); }
	}
}
