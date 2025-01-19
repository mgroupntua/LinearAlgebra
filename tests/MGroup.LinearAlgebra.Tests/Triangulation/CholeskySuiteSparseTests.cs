namespace MGroup.LinearAlgebra.Tests.Triangulation
{
	using MGroup.LinearAlgebra.Implementations.NativeWin64.SuiteSparse;
	using MGroup.LinearAlgebra.Matrices;
	using MGroup.LinearAlgebra.Matrices.Builders;
	using MGroup.LinearAlgebra.Tests;
	using MGroup.LinearAlgebra.Tests.TestData;
	using MGroup.LinearAlgebra.Tests.Utilities;
	using MGroup.LinearAlgebra.Triangulation;
	using MGroup.LinearAlgebra.Vectors;
	using Xunit;

	/// <summary>
	/// Tests for <see cref="CholeskySuiteSparse"/>.
	/// </summary>
	public static class CholeskySuiteSparseTests
	{
		private static readonly MatrixComparer comparer = new MatrixComparer(1E-12);

		[SkippableFact]
		private static void TestRowAddition()
		{
			Skip.IfNot(TestSettings.LibsToTest.Win64SuiteSparse, TestSettings.SkipMessage);

			Matrix original = Matrix.CreateFromArray(SparsePosDef10by10.Matrix);
			Vector rhs = Vector.CreateFromArray(SparsePosDef10by10.Rhs);

			// Start the matrix as diagonal
			var matrixExpected = Matrix.CreateIdentity(original.NumColumns);
			var dok = DokSymmetric.CreateIdentity(SparsePosDef10by10.Order);
			using var factor = new CholeskySuiteSparse(superNodal: false);
			factor.Factorize(dok.BuildSymmetricCscMatrix(true));

			for (int i = 0; i < matrixExpected.NumRows; ++i)
			{
				// Reference solution
				#region minors are not positive definite this way
				//matrixExpected.SetRow(i, newRowVector);
				//matrixExpected.SetColumn(i, newRowVector);
				#endregion
				matrixExpected.SetSubmatrix(0, 0, original.GetSubmatrix(0, i + 1, 0, i + 1)); //this way they are
																							  //Console.WriteLine($"\nOnly dofs [0, {i}]");
																							  //matrixWriter.WriteToConsole(matrixExpected);

				// Update matrix
				Vector newRowVector = matrixExpected.GetRow(i);
				factor.AddRow(i, SparseVector.CreateFromDense(newRowVector));

				// Solve new linear system
				Vector solutionExpected = matrixExpected.FactorCholesky(false).SolveLinearSystem(rhs);
				Vector solutionComputed = factor.SolveLinearSystem(rhs);
				comparer.AssertEqual(solutionExpected, solutionComputed);
			}
		}

		// will probably not work since the matrix will not always be positive definite
		//[SkippableFact]
		private static void TestRowAdditionReverse()
		{
			Skip.IfNot(TestSettings.LibsToTest.Win64SuiteSparse, TestSettings.SkipMessage);

			Matrix original = Matrix.CreateFromArray(SparsePosDef10by10.Matrix);
			Vector rhs = Vector.CreateFromArray(SparsePosDef10by10.Rhs);

			// Start the matrix as diagonal
			var matrixExpected = Matrix.CreateIdentity(original.NumColumns);
			var dok = DokSymmetric.CreateIdentity(SparsePosDef10by10.Order);
			using var factor = new CholeskySuiteSparse(superNodal: false);
			factor.Factorize(dok.BuildSymmetricCscMatrix(true));

			for (int i = 0; i < matrixExpected.NumRows; ++i)
			{
				// Update matrix
				Vector newRowVector = original.GetRow(i);
				matrixExpected.SetSubrow(i, newRowVector);
				matrixExpected.SetSubcolumn(i, newRowVector);
				//Console.WriteLine($"\nOnly dofs [0, {i}]");
				factor.AddRow(i, SparseVector.CreateFromDense(newRowVector));

				// Solve new linear system
				Vector solutionExpected = matrixExpected.FactorCholesky(true).SolveLinearSystem(rhs);
				Vector solutionComputed = factor.SolveLinearSystem(rhs);
				comparer.AssertEqual(solutionExpected, solutionComputed);
			}
		}

		[SkippableFact]
		private static void TestRowDeletion()
		{
			Skip.IfNot(TestSettings.LibsToTest.Win64SuiteSparse, TestSettings.SkipMessage);

			Matrix original = Matrix.CreateFromArray(SparsePosDef10by10.Matrix);
			Vector rhs = Vector.CreateFromArray(SparsePosDef10by10.Rhs);

			// Start the matrix from the original
			var matrixExpected = Matrix.CreateFromArray(SparsePosDef10by10.Matrix);
			var dok = DokSymmetric.CreateEmpty(SparsePosDef10by10.Order);
			for (int j = 0; j < matrixExpected.NumColumns; ++j)
			{
				for (int i = 0; i <= j; ++i)
				{
					if (matrixExpected[i, j] != 0) dok[i, j] = matrixExpected[i, j];
				}
			}

			using var factor = new CholeskySuiteSparse(superNodal: false);
			factor.Factorize(dok.BuildSymmetricCscMatrix(true));

			for (int i = 0; i < matrixExpected.NumRows; ++i)
			{
				// Update matrix
				Vector identityRow = Vector.CreateZero(matrixExpected.NumColumns);
				identityRow[i] = 1.0;
				matrixExpected.SetSubrow(i, identityRow);
				matrixExpected.SetSubcolumn(i, identityRow);
				//Console.WriteLine($"\nOnly dofs [{i + 1}, 10)");
				factor.DeleteRow(i);

				// Solve new linear system
				Vector solutionExpected = matrixExpected.FactorCholesky(false).SolveLinearSystem(rhs);
				Vector solutionComputed = factor.SolveLinearSystem(rhs);
				comparer.AssertEqual(solutionExpected, solutionComputed);
			}
		}

		[SkippableFact]
		private static void TestSystemSolutionSteps()
		{
			Skip.IfNot(TestSettings.LibsToTest.Win64SuiteSparse, TestSettings.SkipMessage);

			int order = SparsePosDef10by10.Order;

			// Build the matrices and right hand sides
			var dense = Matrix.CreateFromArray(SparsePosDef10by10.Matrix);
			//var skyline = SkylineMatrix.CreateFromArrays(order, SparsePositiveDefinite.skylineValues, 
			//    SparsePositiveDefinite.skylineDiagOffsets, false);
			//var dok = DOKSymmetricColMajor.CreateFromSparseMatrix(skyline);
			var dok = DokSymmetric.CreateEmpty(order);
			for (int j = 0; j < order; ++j)
			{
				for (int i = 0; i <= j; ++i)
				{
					if (dense[i, j] != 0) dok[i, j] = dense[i, j];
				}
			}
			Vector b = Vector.CreateFromArray(SparsePosDef10by10.Rhs);
			Matrix B = Matrix.CreateFromArray(SquareInvertible10by10.Matrix);

			// Solve using dense algebra
			CholeskyFull chol = dense.FactorCholesky(false);
			Matrix U = chol.GetFactorU();
			Matrix L = U.Transpose();
			Vector xSolveExpect = chol.SolveLinearSystem(b);
			Matrix XSolveExpect = dense.Invert() * B;
			Vector xBackExpect = U.Invert() * b;
			Matrix XBackExpect = U.Invert() * B;
			Vector xForwardExpect = L.Invert() * b;
			Matrix XForwardExpect = L.Invert() * B;

			// Solve using SuiteSparse
			var (values, rowIndices, colOffsets) = dok.BuildSymmetricCscArrays(true);
			using var factor = new CholeskySuiteSparse(superNodal: true);
			factor.Factorize(order, values.Length, values, rowIndices, colOffsets);
			Vector xSolveComput = factor.SolveLinearSystem(b);
			Matrix XSolveComput = factor.SolveLinearSystems(B);
			Vector xBackComput = factor.BackSubstitution(b);
			Matrix XBackComput = factor.BackSubstitutions(B);
			Vector xForwardComput = factor.ForwardSubstitution(b);
			Matrix XForwardComput = factor.ForwardSubstitutions(B);
			Vector xSolveComput2 = factor.BackSubstitution(factor.ForwardSubstitution(b));

			// Check results
			comparer.AssertEqual(xSolveExpect, xSolveComput);
			comparer.AssertEqual(XSolveExpect, XSolveComput);
			comparer.AssertEqual(xBackExpect, xBackComput);
			comparer.AssertEqual(XBackExpect, XBackComput);
			comparer.AssertEqual(xForwardExpect, xForwardComput);
			comparer.AssertEqual(XForwardExpect, XForwardComput);
		}
	}
}
