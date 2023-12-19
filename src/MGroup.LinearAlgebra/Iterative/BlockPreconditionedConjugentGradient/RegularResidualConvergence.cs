namespace MGroup.LinearAlgebra.Iterative
{
	using System;
	using System.Collections.Generic;
	using System.Diagnostics.CodeAnalysis;
	using System.Security.Cryptography.X509Certificates;
	using System.Text;

	/** Typical residual convergence check.
	 * Nothing fancy with it. Just check if residual^2 < minimal_value
	 */
	public class RegularResidualConvergence : IResidualConvergence
	{
		public double minimal_residual;
		public RegularResidualConvergence(double minimal = 1e-10) { minimal_residual = minimal; }
		public bool isConverged(double rr) => rr < minimal_residual;
	}
}
