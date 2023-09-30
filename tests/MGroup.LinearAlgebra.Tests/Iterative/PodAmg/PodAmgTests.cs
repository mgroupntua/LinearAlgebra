namespace MGroup.LinearAlgebra.Tests.Iterative.PodAmg
{
	using System;
	using System.Collections.Generic;
	using System.Text;

	using MGroup.LinearAlgebra.Iterative;
	using MGroup.LinearAlgebra.Iterative.AlgebraicMultiGrid;
	using MGroup.LinearAlgebra.Iterative.AlgebraicMultiGrid.PodAmg;
	using MGroup.LinearAlgebra.Iterative.AlgebraicMultiGrid.Smoothing;
	using MGroup.LinearAlgebra.Iterative.GaussSeidel;
	using MGroup.LinearAlgebra.Iterative.Termination;
	using MGroup.LinearAlgebra.Matrices;
	using MGroup.LinearAlgebra.Tests.Utilities;
	using MGroup.LinearAlgebra.Vectors;

	using Xunit;

	public static class PodAmgTests
	{
		[Fact]
		public static void TestPodAMG()
		{
			var matrix = Matrix.CreateFromArray(DataSet1.Matrix);
			var csr = CsrMatrix.CreateFromDense(matrix);
			var rhs = Vector.CreateFromArray(DataSet1.Rhs);
			var solutionExpected = Vector.CreateFromArray(DataSet1.Solution);
			var samples = Matrix.CreateFromArray(DataSet1.Samples);
			int numPrincipalComponents = DataSet1.PrincipalComponents.GetLength(1);
			int numPodAmgCyclesExpected = DataSet1.NumPodAmgCycles;

			var solverBuilder = new PodAmgAlgorithm.Builder();
			solverBuilder.MaxIterationsProvider = new FixedMaxIterationsProvider(30000);
			solverBuilder.ConvergenceTolerance = 1E-5;
			solverBuilder.ConvergenceCriterion = new SolutionNeverConvergesCriterion();
			solverBuilder.SmootherBuilder = new GaussSeidelSmoother.Builder(
				new GaussSeidelIterationCsrSerial.Builder(), GaussSeidelSweepDirection.Symmetric, numIterations: 1);

			using PodAmgAlgorithm solver = solverBuilder.Create(samples, numPrincipalComponents);
			var solutionComputed = Vector.CreateZero(rhs.Length);
			solver.Initialize(csr);
			IterativeStatistics stats = solver.Solve(rhs, solutionComputed);

			var comparer = new MatrixComparer(1E-6);
			comparer.AssertEqual(solutionExpected, solutionComputed);
			Assert.InRange(stats.NumIterationsRequired, 0, numPodAmgCyclesExpected);
		}

	}
}
