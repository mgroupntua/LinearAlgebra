using System.Text;

using MGroup.LinearAlgebra.Iterative;
using MGroup.LinearAlgebra.Iterative.GaussSeidel;
using MGroup.LinearAlgebra.Iterative.Termination.Convergence;
using MGroup.LinearAlgebra.Iterative.Termination.Iterations;
using MGroup.LinearAlgebra.Matrices;
using MGroup.LinearAlgebra.Output;
using MGroup.LinearAlgebra.Tests.TestData;
using MGroup.LinearAlgebra.Tests.Utilities;
using MGroup.LinearAlgebra.Vectors;
using Xunit;

namespace MGroup.LinearAlgebra.Tests.Iterative
{
	/// <summary>
	/// Tests for <see cref="GaussSeidel"/>.
	/// Authors: Gerasimos Sotiropoulos 
	/// </summary>
	public static class GaussSeidelTests
	{
		[Theory]
		[InlineData(true, 10, 1E-9, 1E-5, true)]
		[InlineData(false, 10, 1E-8, 1E-5, true)]
		[InlineData(true, 10, 1E-9, 1E-5, false)]
		[InlineData(false, 10, 1E-8, 1E-5, false)]
		private static void TestSparseSystem(bool forwardGaussSeidel, int numIterations, double gsConvergenceTolerance, double entrywiseTolerance, bool csrFormat)
		{
			IMatrixView matrix;
			IGaussSeidelIteration gsIteration;
			if (csrFormat)
			{
				matrix = CsrMatrix.CreateFromArrays(SparsePosDef10by10.Order, SparsePosDef10by10.Order, SparsePosDef10by10.CsrValues, 
					SparsePosDef10by10.CsrColIndices, SparsePosDef10by10.CsrRowOffsets, true);
				gsIteration = new GaussSeidelIterationCsrSerial();
			}
			else
			{
				matrix = Matrix.CreateFromArray(SparsePosDef10by10.Matrix);
				gsIteration = new GaussSeidelIterationGeneral();
			}

			var b = Vector.CreateFromArray(SparsePosDef10by10.Rhs);
			var xExpected = Vector.CreateFromArray(SparsePosDef10by10.Lhs);

			var builder = new GaussSeidelAlgorithm.Builder();
			//builder.ConvergenceCriterion = new SolutionNeverConvergesCriterion(); // We would use this, but to test we need to track the convergence rate.
			builder.ConvergenceCriterion = new AbsoluteSolutionConvergenceCriterion();
			builder.ConvergenceTolerance = 0.0;
			builder.MaxIterationsProvider = new FixedMaxIterationsProvider(numIterations);
			builder.ForwardGaussSeidel = forwardGaussSeidel;
			var gs = builder.Build(gsIteration);
			var xComputed = Vector.CreateZero(b.Length);

			IterativeStatistics stats = gs.Solve(matrix, b, xComputed);

			Assert.Equal(numIterations, stats.NumIterationsRequired);
			Assert.InRange(stats.ConvergenceCriterion.value, 0.0, gsConvergenceTolerance);

			var comparer = new MatrixComparer(entrywiseTolerance);
			comparer.AssertEqual(xExpected, xComputed);


		}
	}
}
