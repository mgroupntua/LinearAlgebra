using MGroup.LinearAlgebra.Iterative;
using MGroup.LinearAlgebra.Iterative.Termination;

namespace MGroup.LinearAlgebra.Distributed.IterativeMethods.PCG
{
	public abstract class PcgBuilderBase
	{
		/// <summary>
		/// Specifies how to calculate the maximum iterations that the PCG algorithm will run for.
		/// </summary>
		public IMaxIterationsProvider MaxIterationsProvider { get; set; } = new PercentageMaxIterationsProvider(1.0);

		/// <summary>
		/// Specifies how the PCG algorithm will check that convergence has been reached.
		/// </summary>
		public IPcgResidualConvergence Convergence { get; set; } = new RegularPcgConvergence();

		/// <summary>
		/// Specifies how often the residual vector will be corrected by an exact (but costly) calculation.
		/// </summary>
		public IPcgResidualUpdater ResidualUpdater { get; set; } = new RegularPcgResidualUpdater();

		/// <summary>
		/// The PCG algorithm will converge when sqrt(r*inv(M)*r) / sqrt(r0*inv(M)*r0) &lt;= <paramref name="ResidualTolerance"/>,
		/// where M is the preconditioner, r = A*x is the current residual vector and r0 = A*x0 the initial residual vector.
		/// </summary>
		public double ResidualTolerance { get; set; } = 1E-10;

		/// <summary>
		/// If the iterative algorithm does not converge to the specified tolerance <see cref="ResidualTolerance"/>, within
		/// the number of iterations specified by <see cref="MaxIterationsProvider"/>, then the algorithm will terminate.
		/// If <see cref="ThrowExceptionIfNotConvergence"/> == true, then an exception will be thrown. Else execution will be 
		/// resumed by the client and the returned <see cref="IterativeStatistics"/> object will report that convergence was not
		/// reached.
		/// </summary>
		public bool ThrowExceptionIfNotConvergence { get; set; } = true; 
		//TODO: perhaps it would be better to allow the clients to specify how convergence is defined. E.g. only iterations,
		//		only residual, both, some other criterion etc.
	}
}
