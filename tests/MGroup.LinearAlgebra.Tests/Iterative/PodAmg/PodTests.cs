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
		[InlineData(false)]
		public static void TestPod(bool samplesAreMoreThanProblemDimension)
		{
			Matrix samples, principalComponentsExpected;
			if (samplesAreMoreThanProblemDimension)
			{
				throw new NotImplementedException();
				samples = Matrix.CreateFromArray(DataSet2.Samples);
				principalComponentsExpected = Matrix.CreateFromArray(DataSet2.PrincipalComponents);
			}
			else
			{
				samples = Matrix.CreateFromArray(DataSet1.Samples);
				principalComponentsExpected = Matrix.CreateFromArray(DataSet1.PrincipalComponents);
			}

			int numSamples = samples.NumColumns;
			int numPrincipalComponents = principalComponentsExpected.NumColumns;

			var pod = new ProperOrthogonalDecomposition();
			Matrix principalComponentsComputed = pod.CalculatePrincipalComponents(numSamples, samples, numPrincipalComponents);

			var comparer = new MatrixComparer(1E-12);
			comparer.AssertEqual(principalComponentsExpected, principalComponentsComputed);
		}
	}
}
