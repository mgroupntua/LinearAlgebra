using System;
using System.Collections.Generic;
using System.Text;
using Xunit;
using MGroup.LinearAlgebra.Vectors;
using MGroup.LinearAlgebra.Matrices;
using MGroup.LinearAlgebra.Tests.TestData;
using MGroup.LinearAlgebra.Eigensystems;
using MGroup.LinearAlgebra.Tests.Utilities;

namespace MGroup.LinearAlgebra.Tests.Eigensystems
{
	public static class SymmetricEigensystemFullTests
	{
		private static readonly MatrixComparer comparer = new MatrixComparer(1E-10);

		[Fact]
		private static void TestFullEigendecomposition()
		{
			double tol = 1E-10;
			var A = Matrix.CreateFromArray(SymmPosDef10by10.Matrix);

			(Vector lambda, Matrix X) = SpectralUtilities.SortSingularValues(
				Vector.CreateFromArray(SymmPosDef10by10.Eigenvalues),
				Matrix.CreateFromArray(SymmPosDef10by10.Eigenvectors),
				descending: true);
			var decompExpected = new SpectralDecomposition(lambda, X, tol);

			var eig = SymmetricEigensystemFull.Create(A.NumColumns, A.Copy().RawData, true);
			var decompComputed = new SpectralDecomposition(eig.EigenvaluesReal, eig.EigenvectorsRight, tol);

			Assert.True(decompExpected.IsEquivalent(decompComputed));
			Assert.True(decompComputed.CanRecomposeOriginalMatrix(A));

			// Very low chance that a libray outputs eigenvectors normalized with a constistent sign
			//Assert.True(decompExpected.IsIdentical(decompComputed)); 
		}

		[Theory]
		[MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
		private static void TestEigenvaluesAndEigenvectors(LinearAlgebraProviderChoice providers)
		{
			TestSettings.RunMultiproviderTest(providers, delegate ()
			{
				var A = Matrix.CreateFromArray(SymmPosDef10by10.Matrix);
				var eigenvaluesExpected = Vector.CreateFromArray(SymmPosDef10by10.Eigenvalues);
				var eigenvectorsExpected = Matrix.CreateFromArray(SymmPosDef10by10.Eigenvectors);
				var eigensystem = SymmetricEigensystemFull.Create(A.NumRows, A.RawData, true);

				// Check
				comparer.AssertEqual(eigenvaluesExpected, eigensystem.EigenvaluesReal);
				comparer.AssertEqual(eigenvectorsExpected, eigensystem.EigenvectorsRight);
				Assert.True(eigensystem.EigenvaluesImaginary == null);
				Assert.True(eigensystem.EigenvectorsLeft == null);
			});
		}

		[Theory]
		[MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
		private static void TestEigenvaluesOnly(LinearAlgebraProviderChoice providers)
		{
			TestSettings.RunMultiproviderTest(providers, delegate ()
			{
				var A = Matrix.CreateFromArray(SymmPosDef10by10.Matrix);
				var eigenvaluesExpected = Vector.CreateFromArray(SymmPosDef10by10.Eigenvalues);
				var eigenvectorsExpected = Matrix.CreateFromArray(SymmPosDef10by10.Eigenvectors);
				var eigensystem = SymmetricEigensystemFull.Create(A.NumRows, A.RawData, false);

				// Check
				comparer.AssertEqual(eigenvaluesExpected, eigensystem.EigenvaluesReal);
				Assert.True(eigensystem.EigenvectorsRight == null);
				Assert.True(eigensystem.EigenvaluesImaginary == null);
				Assert.True(eigensystem.EigenvectorsLeft == null);
			});
		}

		[Theory]
		[MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
		private static void TestEigensystemCalledFromFullMatrix(LinearAlgebraProviderChoice providers)
		{
			TestSettings.RunMultiproviderTest(providers, delegate ()
			{
				var A = Matrix.CreateFromArray(SymmPosDef10by10.Matrix);
				var eigenvaluesExpected = Vector.CreateFromArray(SymmPosDef10by10.Eigenvalues);
				var eigenvectorsExpected = Matrix.CreateFromArray(SymmPosDef10by10.Eigenvectors);
				(Vector eigenvalues, Matrix eigenvectors) = A.CalcEigensystemSymmetric();

				// Check
				comparer.AssertEqual(eigenvaluesExpected, eigenvalues);
				comparer.AssertEqual(eigenvectorsExpected, eigenvectors);
			});
		}

		[Theory]
		[MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
		private static void TestEigensystemCalledFromPackedMatrix(LinearAlgebraProviderChoice providers)
		{
			TestSettings.RunMultiproviderTest(providers, delegate ()
			{
				var A = SymmetricMatrix.CreateFromArray(SymmPosDef10by10.Matrix);
				var eigenvaluesExpected = Vector.CreateFromArray(SymmPosDef10by10.Eigenvalues);
				var eigenvectorsExpected = Matrix.CreateFromArray(SymmPosDef10by10.Eigenvectors);
				(Vector eigenvalues, Matrix eigenvectors) = A.CalcEigensystem();

				// Check
				comparer.AssertEqual(eigenvaluesExpected, eigenvalues);
				comparer.AssertEqual(eigenvectorsExpected, eigenvectors);
			});
		}
	}
}
