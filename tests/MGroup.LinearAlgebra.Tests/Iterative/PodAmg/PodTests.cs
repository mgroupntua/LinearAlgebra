namespace MGroup.LinearAlgebra.Tests.Iterative.PodAmg
{
	using System;
	using System.Collections.Generic;
	using System.Text;

	using MGroup.LinearAlgebra.Iterative.AlgebraicMultiGrid.PodAmg;
	using MGroup.LinearAlgebra.Matrices;
	using MGroup.LinearAlgebra.Tests.Utilities;

	using Xunit;

	public static class PodTests
	{
		[Theory]
		[InlineData(false, false)]
		[InlineData(false, true)]
		public static void TestPod(bool samplesAreMoreThanProblemDimension, bool keepOnlyNonZeroEigenvalues)
		{
			Matrix samples, principalComponentsExpected;
			double tol;
			if (samplesAreMoreThanProblemDimension)
			{
				throw new NotImplementedException();
				samples = Matrix.CreateFromArray(DataSet3.Samples);
				principalComponentsExpected = Matrix.CreateFromArray(DataSet3.PrincipalComponents);
				tol = 0;
			}
			else
			{
				if (keepOnlyNonZeroEigenvalues)
				{
					samples = Matrix.CreateFromArray(DataSet2.Samples);
					principalComponentsExpected = Matrix.CreateFromArray(DataSet2.PrincipalComponents);
					tol = 1E-11;
				}
				else
				{
					samples = Matrix.CreateFromArray(DataSet1.Samples);
					principalComponentsExpected = Matrix.CreateFromArray(DataSet1.PrincipalComponents);
					tol = 1E-12;
				}
			}

			int numSamples = samples.NumColumns;
			int numPrincipalComponents = principalComponentsExpected.NumColumns;

			var pod = new ProperOrthogonalDecomposition(keepOnlyNonZeroEigenvalues);
			Matrix principalComponentsComputed = pod.CalculatePrincipalComponents(numSamples, samples, numPrincipalComponents);

			var comparer = new MatrixComparer(tol);
			comparer.AssertEqual(principalComponentsExpected, principalComponentsComputed);
		}
	}
}
