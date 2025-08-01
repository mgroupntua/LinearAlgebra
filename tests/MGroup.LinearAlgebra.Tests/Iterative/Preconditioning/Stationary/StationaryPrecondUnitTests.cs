namespace MGroup.LinearAlgebra.Tests.Iterative.Preconditioning.Stationary
{
	using System;
	using System.Collections;
	using System.Collections.Generic;

	using MGroup.LinearAlgebra.Iterative.Preconditioning;
	using MGroup.LinearAlgebra.Iterative.Preconditioning.Stationary;
	using MGroup.LinearAlgebra.Matrices;
	using MGroup.LinearAlgebra.Reduction;
	using MGroup.LinearAlgebra.Tests.TestData;
	using MGroup.LinearAlgebra.Tests.Utilities;
	using MGroup.LinearAlgebra.Vectors;

	using Xunit;

	public class StationaryPrecondUnitTests
	{
		[Fact]
		public static void TestGaussSeidelBackPreconditioner()
		{
			var preconditioner = new GaussSeidelPreconditionerCsr(forwardDirection: false);
			RunTwoApplications(preconditioner, (decomp, rhs) =>
			{
				// Solve (D+U) * x = y
				var UplusD = decomp.GetCombination("U+D", 0, 1, 1);
				return UplusD.SolveBackSubstitution(rhs);
			});
		}

		[Fact]
		public static void TestGaussSeidelForwardPreconditioner()
		{
			var preconditioner = new GaussSeidelPreconditionerCsr(forwardDirection: true);
			RunTwoApplications(preconditioner, (decomp, rhs) =>
			{
				// Solve (D+L) * x = y
				var LplusD = decomp.GetCombination("L+D", 1, 1, 0);
				return LplusD.SolveForwardSubstitution(rhs);
			});
		}

		[Fact]
		public static void TestJacobiPreconditioner()
		{
			var preconditioner = new JacobiPreconditioner(1E-10);
			RunTwoApplications(preconditioner, (decomp, rhs) =>
			{
				// Solve D * x = y
				var D = decomp.GetD();
				return D.SolveDiagonal(rhs);
			});
		}

		[Fact]
		public static void TestJacobiPreconditioner1Application()
		{
			var comparer = new MatrixComparer(1E-10);
			var matrix = Matrix.CreateFromArray(SquareSingular10by10.Matrix);
			var preconditioner = new JacobiPreconditioner();
			preconditioner.UpdateMatrix(matrix, true);

			var b = Vector.CreateWithValue(10, 1.0);
			var xExpected = Vector.CreateFromArray(new double[]
			{
				1, 3.000300030003000, 0.25, 0.2, 4, 0.111111111111111, 0.125, 0.142857142857143, 1, 0.5
			});
			var xComputed = Vector.CreateZero(10);

			preconditioner.SolveLinearSystem(b, xComputed);
			comparer.AssertEqual(xExpected, xComputed);
		}

		[Fact]
		public static void TestSorBackPreconditioner()
		{
			var omega = 1.2;
			var preconditioner = new SorPreconditionerCsr(omega, forwardDirection: false);
			RunTwoApplications(preconditioner, (decomp, rhs) =>
			{
				// Solve 1/omega*(D+omega*U) * x = y
				var DplusWU = decomp.GetCombination("D+ωU", 0, 1, omega);
				return DplusWU.SolveBackSubstitution(rhs).Scale(omega);
			});
		}

		[Fact]
		public static void TestSorForwardPreconditioner()
		{
			var omega = 1.2;
			var preconditioner = new SorPreconditionerCsr(omega, forwardDirection: true);
			RunTwoApplications(preconditioner, (decomp, rhs) =>
			{
				// Solve 1/omega*(D+omega*L) * x = y
				var DplusWL = decomp.GetCombination("D+ωL", omega, 1, 0);
				return DplusWL.SolveForwardSubstitution(rhs).Scale(omega);
			});
		}

		[Fact]
		public static void TestSsorPreconditioner()
		{
			var omega = 1.2;
			var preconditioner = new SsorPreconditionerCsr(omega);
			RunTwoApplications(preconditioner, (decomp, rhs) =>
			{
				// Solve 1/(omega*(2-omega)*(D+omega*L)*inv(D)*(D+omega*U) * x = y
				var D = decomp.GetD();
				var DplusWL = decomp.GetCombination("D+ωL", omega, 1, 0);
				var DplusWU = decomp.GetCombination("D+ωU", 0, 1, omega);

				var x0 = rhs.Scale(omega * (2 - omega));
				var x1 = DplusWL.SolveForwardSubstitution(x0);
				var x2 = D * x1;
				return DplusWU.SolveBackSubstitution(x2);
			});
		}

		private static void RunTwoApplications(IPreconditioner preconditioner,
			Func<StationaryAlgorithmDecomposition, Vector, Vector> preconditionWithMatrixForm)
		{
			// Setup comparison code
			var entrywiseTolerance = 1E-15;
			var comparer = new MatrixComparer(entrywiseTolerance);

			// Initialize matrices and vectors
			var csrMatrix = CsrMatrix.CreateFromArrays(SparsePosDef10by10.Order, SparsePosDef10by10.Order,
					SparsePosDef10by10.CsrValues, SparsePosDef10by10.CsrColIndices, SparsePosDef10by10.CsrRowOffsets, true);
			var y1 = Vector.CreateFromArray(SparsePosDef10by10.Rhs);
			var max = y1.MaxAbsolute();
			var min = y1.MinAbsolute();
			var y2 = y1.Scale(0.1);
			y2.DoToAllEntriesIntoThis(x => x + 0.01 * (max - min) / max);
			var xExpected = Vector.CreateZero(y1.Length);
			var xComputed = Vector.CreateZero(y1.Length);

			// Prepare preconditioner
			var decomp = new StationaryAlgorithmDecomposition(SparseMatrix.CreateFromMatrix(csrMatrix));
			preconditioner.UpdateMatrix(csrMatrix, true);

			// First try
			xExpected = preconditionWithMatrixForm(decomp, y1);
			preconditioner.SolveLinearSystem(y1, xComputed);
			comparer.AssertEqual(xExpected, xComputed);

			// Second try
			xExpected = preconditionWithMatrixForm(decomp, y2);
			preconditioner.SolveLinearSystem(y2, xComputed);
			comparer.AssertEqual(xExpected, xComputed);
		}
	}
}
