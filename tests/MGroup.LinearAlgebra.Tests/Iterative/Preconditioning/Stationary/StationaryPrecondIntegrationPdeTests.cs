namespace MGroup.LinearAlgebra.Tests.Iterative.Preconditioning.Stationary
{
	using MGroup.LinearAlgebra.Iterative.GeneralizedMinimalResidual;
	using MGroup.LinearAlgebra.Iterative.PreconditionedConjugateGradient;
	using MGroup.LinearAlgebra.Iterative.Preconditioning;
	using MGroup.LinearAlgebra.Iterative.Preconditioning.Stationary;
	using MGroup.LinearAlgebra.Iterative.Termination.Iterations;
	using MGroup.LinearAlgebra.Matrices;
	using MGroup.LinearAlgebra.Matrices.Builders;
	using MGroup.LinearAlgebra.Tests.TestData;
	using MGroup.LinearAlgebra.Tests.TestData.SparseLinearSystems;
	using MGroup.LinearAlgebra.Tests.Utilities;
	using MGroup.LinearAlgebra.Vectors;

	using Xunit;

	public class StationaryPrecondIntegrationPdeTests
	{
		[Fact]
		private static void TestJacobiGmresFdm2D()
		{
			SolveProblemFd2D(new JacobiPreconditioner(1E-10), 1E-7, 276);
		}

		[Fact]
		private static void TestJacobiPcgFem2D()
		{
			SolveProblemFem2D(new JacobiPreconditioner(1E-10), 1E-7, 86);
		}

		[Fact]
		private static void TestJacobiPcgFem3D()
		{
			SolveProblemFem3D(new JacobiPreconditioner(1E-10), 1E-7, 138);
		}

		[Theory]
		[InlineData(true, 218)]
		[InlineData(false, 223)]
		private static void TestGaussSeidelGmresFdm2D(bool forwardDirection, int numIterationsExpected)
		{
			SolveProblemFd2D(new GaussSeidelPreconditionerCsr(forwardDirection, numApplications: 1), 1E-7, numIterationsExpected);
		}

		[Theory]
		[InlineData(true, 205)]
		[InlineData(false, 206)]
		private static void TestSorGmresFdm2D(bool forwardDirection, int numIterationsExpected)
		{
			SolveProblemFd2D(new SorPreconditionerCsr(1.8, forwardDirection, numApplications: 1), 1E-7, numIterationsExpected);
		}

		[Fact]
		private static void TestSsorGmresFdm2D()
		{
			SolveProblemFd2D(new SsorPreconditionerCsr(1.9, numApplications: 1), 1E-7, 37);
		}

		[Fact]
		private static void TestSsorPcgFem2D()
		{
			SolveProblemFem2D(new SsorPreconditionerCsr(1.2, numApplications: 1), 1E-7, 43);
		}

		[Fact]
		private static void TestSsorPcgFem3D()
		{
			SolveProblemFem3D(new SsorPreconditionerCsr(1.2, numApplications: 1), 1E-7, 104);
		}

		private static void SolveProblemFd2D(IPreconditioner preconditioner, double residualTol, int numIterations)
		{
			var model = new FdmPoisson2D(101);
			(SparseMatrix A, Vector xExpected, Vector b) = model.CreateLinearSystem();
			SolveWithGmres(A, xExpected, b, preconditioner, residualTol, numIterations);
		}

		private static void SolveProblemFem2D(IPreconditioner preconditioner, double residualTol, int numIterations)
		{
			var mesh = new CartesianMesh2D(25, 5, 10.0, 2.0);
			var model = new FemCantilever2D(mesh, 1.0);
			model.ElasticityModulus = 2E7;
			(SparseMatrix A, Vector xExpected, Vector b) = model.CreateLinearSystem();
			SolveWithPcg(A, xExpected, b, preconditioner, residualTol, numIterations);
		}

		private static void SolveProblemFem3D(IPreconditioner preconditioner, double residualTol, int numIterations)
		{
			var mesh = new CartesianMesh3D(25, 5, 5, 10.0, 1.0, 2.0);
			var model = new FemCantilever3D(mesh);
			model.ElasticityModulus = 2E7;
			(SparseMatrix A, Vector xExpected, Vector b) = model.CreateLinearSystem();
			SolveWithPcg(A, xExpected, b, preconditioner, residualTol, numIterations);
		}

		private static void SolveWithGmres(SparseMatrix A, Vector xExpected, Vector b,
			IPreconditioner preconditioner, double residualTol, int numIterations)
		{
			var csrA = DokRowMajor.CreateFromSparsePattern(A.NumRows, A.NumColumns, A.EnumerateNonZeros()).BuildCsrMatrix(true);
			var xComputed = Vector.CreateZero(xExpected.Length);

			var builder = new GmresAlgorithm.Builder();
			builder.RelativeTolerance = residualTol;
			builder.AbsoluteTolerance = residualTol * 100;
			builder.MaximumIterations = 10 * A.NumRows;
			var gmres = builder.Build();
			preconditioner.UpdateMatrix(csrA, true);
			var stats = gmres.Solve(csrA, preconditioner, b, xComputed, true);

			var comparer = new MatrixComparer(1E-5);
			comparer.AssertEqual(xExpected, xComputed);
			Assert.InRange(stats.NumIterationsRequired, 0, numIterations);
		}

		private static void SolveWithPcg(SparseMatrix A, Vector xExpected, Vector b,
			IPreconditioner preconditioner, double residualTol, int numIterations)
		{
			var csrA = DokRowMajor.CreateFromSparsePattern(A.NumRows, A.NumColumns, A.EnumerateNonZeros()).BuildCsrMatrix(true);
			var xComputed = Vector.CreateZero(xExpected.Length);

			var builder = new PcgAlgorithm.Factory();
			builder.ResidualTolerance = residualTol;
			builder.MaxIterationsProvider = new PercentageMaxIterationsProvider(1.0);
			var pcg = builder.Build();
			preconditioner.UpdateMatrix(csrA, true);
			var stats = pcg.Solve(csrA, preconditioner, b, xComputed, true);

			var comparer = new MatrixComparer(1E-5);
			comparer.AssertEqual(xExpected, xComputed);
			Assert.InRange(stats.NumIterationsRequired, 0, numIterations);
		}
	}
}
