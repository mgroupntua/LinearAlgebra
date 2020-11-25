using System;
using System.Collections.Generic;
using System.Text;
using Xunit;
using MGroup.LinearAlgebra.Vectors;
using MGroup.LinearAlgebra.Matrices;
using MGroup.LinearAlgebra.Tests.TestData;
using MGroup.LinearAlgebra.Eigensystems;
using MGroup.LinearAlgebra.Tests.Utilities;

//TODO: Test all 4 combos (only values, all three, values and left or right
//TODO: Test for matrices that have all real eigenvalues and for matrices that have complex
//TODO: Also test for Matrix.CalcEigenSystemNonSymmetric()
namespace MGroup.LinearAlgebra.Tests.Eigensystems
{
	public static class NonSymmetricEigensystemFullTests
	{
		private static readonly MatrixComparer comparer = new MatrixComparer(1E-10);

		[Theory]
		[MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
		private static void TestEigenvaluesOnly(LinearAlgebraProviderChoice providers)
		{
			TestSettings.RunMultiproviderTest(providers, delegate ()
			{
				var A = Matrix.CreateFromArray(SquareInvertible10by10.Matrix);
				var eigenvaluesRealExpected = Vector.CreateFromArray(SquareInvertible10by10.EigenvaluesRealPart);
				var eigenvaluesImaginaryExpected = Vector.CreateFromArray(SquareInvertible10by10.EigenvaluesImaginaryPart);
				var eigensystem = NonSymmetricEigensystemFull.Create(A.NumRows, A.RawData, false, false);

				// Check
				comparer.AssertEqual(eigenvaluesRealExpected, eigensystem.EigenvaluesReal);
				comparer.AssertEqual(eigenvaluesImaginaryExpected, eigensystem.EigenvaluesImaginary);
				Assert.True(eigensystem.EigenvectorsLeft == null);
				Assert.True(eigensystem.EigenvectorsRight == null);
			});
		}

		[Theory]
		[MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
		private static void TestEigenvaluesAndRightEigenvectors(LinearAlgebraProviderChoice providers)
		{
			TestSettings.RunMultiproviderTest(providers, delegate ()
			{
				var A = Matrix.CreateFromArray(SquareInvertible10by10.Matrix);
				var eigenvaluesRealExpected = Vector.CreateFromArray(SquareInvertible10by10.EigenvaluesRealPart);
				var eigenvaluesImaginaryExpected = Vector.CreateFromArray(SquareInvertible10by10.EigenvaluesImaginaryPart);
				var eigenvectorsRightExpected = Matrix.CreateFromArray(SquareInvertible10by10.EigenvectorsRight);
				var eigensystem = NonSymmetricEigensystemFull.Create(A.NumRows, A.RawData, false, true);

				// Check
				comparer.AssertEqual(eigenvaluesRealExpected, eigensystem.EigenvaluesReal);
				comparer.AssertEqual(eigenvaluesImaginaryExpected, eigensystem.EigenvaluesImaginary);
				comparer.AssertEqual(eigenvectorsRightExpected, eigensystem.EigenvectorsRight);
				Assert.True(eigensystem.EigenvectorsLeft == null);
			});
		}

		[Theory]
		[MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
		private static void TestEigenvaluesAndLeftEigenvectors(LinearAlgebraProviderChoice providers)
		{
			TestSettings.RunMultiproviderTest(providers, delegate ()
			{
				var A = Matrix.CreateFromArray(SquareInvertible10by10.Matrix);
				var eigenvaluesRealExpected = Vector.CreateFromArray(SquareInvertible10by10.EigenvaluesRealPart);
				var eigenvaluesImaginaryExpected = Vector.CreateFromArray(SquareInvertible10by10.EigenvaluesImaginaryPart);
				var eigenvectorsLeftExpected = Matrix.CreateFromArray(SquareInvertible10by10.EigenvectorsLeft);
				var eigensystem = NonSymmetricEigensystemFull.Create(A.NumRows, A.RawData, true, false);

				// Check
				comparer.AssertEqual(eigenvaluesRealExpected, eigensystem.EigenvaluesReal);
				comparer.AssertEqual(eigenvaluesImaginaryExpected, eigensystem.EigenvaluesImaginary);
				comparer.AssertEqual(eigenvectorsLeftExpected, eigensystem.EigenvectorsLeft);
				Assert.True(eigensystem.EigenvectorsRight == null);
			});
		}

		[Theory]
		[MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
		private static void TestEigenvaluesAndAllEigenvectors(LinearAlgebraProviderChoice providers)
		{
			TestSettings.RunMultiproviderTest(providers, delegate ()
			{
				var A = Matrix.CreateFromArray(SquareInvertible10by10.Matrix);
				var eigenvaluesRealExpected = Vector.CreateFromArray(SquareInvertible10by10.EigenvaluesRealPart);
				var eigenvaluesImaginaryExpected = Vector.CreateFromArray(SquareInvertible10by10.EigenvaluesImaginaryPart);
				var eigenvectorsRightExpected = Matrix.CreateFromArray(SquareInvertible10by10.EigenvectorsRight);
				var eigenvectorsLeftExpected = Matrix.CreateFromArray(SquareInvertible10by10.EigenvectorsLeft);
				var eigensystem = NonSymmetricEigensystemFull.Create(A.NumRows, A.RawData, true, true);

				// Check
				comparer.AssertEqual(eigenvaluesRealExpected, eigensystem.EigenvaluesReal);
				comparer.AssertEqual(eigenvaluesImaginaryExpected, eigensystem.EigenvaluesImaginary);
				comparer.AssertEqual(eigenvectorsRightExpected, eigensystem.EigenvectorsRight);
				comparer.AssertEqual(eigenvectorsLeftExpected, eigensystem.EigenvectorsLeft);
			});
		}

		[Theory]
		[MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
		private static void TestEigensystemCalledFromFullMatrix(LinearAlgebraProviderChoice providers)
		{
			TestSettings.RunMultiproviderTest(providers, delegate ()
			{
				var A = Matrix.CreateFromArray(SquareInvertible10by10.Matrix);
				var eigenvaluesRealExpected = Vector.CreateFromArray(SquareInvertible10by10.EigenvaluesRealPart);
				var eigenvaluesImaginaryExpected = Vector.CreateFromArray(SquareInvertible10by10.EigenvaluesImaginaryPart);
				var eigenvectorsRightExpected = Matrix.CreateFromArray(SquareInvertible10by10.EigenvectorsRight);
				var eigenvectorsLeftExpected = Matrix.CreateFromArray(SquareInvertible10by10.EigenvectorsLeft);
				(Vector eigenvaluesReal, Vector eigenvaluesImaginary, Matrix eigenvectorsRight, Matrix eigenvectorsLeft) =
					A.CalcEigensystemNonSymmetric();

				// Check
				comparer.AssertEqual(eigenvaluesRealExpected, eigenvaluesReal);
				comparer.AssertEqual(eigenvaluesImaginaryExpected, eigenvaluesImaginary);
				comparer.AssertEqual(eigenvectorsRightExpected, eigenvectorsRight);
				comparer.AssertEqual(eigenvectorsLeftExpected, eigenvectorsLeft);
			});
		}
	}
}
