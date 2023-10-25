namespace MGroup.LinearAlgebra.Tests.Eigensystems
{
	using System;
	using System.Collections;
	using System.Collections.Generic;
	using System.Linq;
	using System.Reflection;
	using System.Text;
	using MGroup.LinearAlgebra.Eigensystems;
	using MGroup.LinearAlgebra.Matrices;
	using MGroup.LinearAlgebra.Tests.TestData;

	using MGroup.LinearAlgebra.Tests.Utilities;
	using MGroup.LinearAlgebra.Vectors;

	using Xunit;

	public class SingularValueDecompositionTests
	{
		[Fact]
		private static void TestSingluarValuesAndEigenvectors()
		{
			double tol = 1E-10;
			var A = Matrix.CreateFromArray(SymmPosDef10by10.Matrix);

			(Vector lambda, Matrix X) = SpectralUtilities.SortSingularValues(
				Vector.CreateFromArray(SymmPosDef10by10.Eigenvalues),
				Matrix.CreateFromArray(SymmPosDef10by10.Eigenvectors), 
				descending:true);
			var decompExpected = new SpectralDecomposition(lambda, X, tol);

			var svd = SingularValueDecomposition.Calculate(A);
			var decompComputed = new SpectralDecomposition(svd.SingularValues, svd.SingularVectors, tol);

			Assert.True(decompExpected.IsEquivalent(decompComputed));
			Assert.True(decompComputed.CanRecomposeOriginalMatrix(A));

			// Very low chance that a libray outputs eigenvectors normalized with a constistent sign
			//Assert.True(decompExpected.IsIdentical(decompComputed)); 
		}

	}
}
