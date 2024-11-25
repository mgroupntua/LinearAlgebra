namespace MGroup.LinearAlgebra.Tests.Reordering
{
	using System;
	using System.Collections.Generic;
	using System.Linq;
	using System.Text;
	using System.Threading.Tasks;

	using MGroup.LinearAlgebra.Implementations;
	using MGroup.LinearAlgebra.Implementations.Managed;
	using MGroup.LinearAlgebra.Matrices;
	using MGroup.LinearAlgebra.Reordering;
	using MGroup.LinearAlgebra.Tests.TestData;
	using MGroup.LinearAlgebra.Tests.Utilities;

	using Xunit;

	/// <summary>
	/// Tests for <see cref="AmdSymmetricOrdering"/>.
	/// </summary>
	public class AmdSymmetricOrderingTests
	{
		private static readonly MatrixComparer comparer = new MatrixComparer(1E-13);

		[Theory]
		[MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
		public static void TestFindPermutationGivenPattern(IImplementationProvider provider)
		{
			TestSettings.RunMultiproviderTest(provider, delegate ()
			{
				var pattern = SparsityPatternSymmetric.CreateFromDense(Matrix.CreateFromArray(SparsePosDef10by10.Matrix));
				var orderingAlg = new AmdSymmetricOrdering();
				(int[] permutation, bool oldToNew) = orderingAlg.FindPermutation(pattern);
				Assert.True(!oldToNew);
				comparer.AssertEqual(SparsePosDef10by10.MatlabPermutationAMD, permutation);
			});
		}

		[Theory]
		[MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
		public static void TestFindPermutationGivenRawArrays(IImplementationProvider provider)
		{
			TestSettings.RunMultiproviderTest(provider, delegate ()
			{
				int order = SparsePosDef10by10.Order;
				int[] cscColOffsets = SparsePosDef10by10.SymmetricCscColOffsets;
				int[] cscRowIndices = SparsePosDef10by10.SymmetricCscRowIndices;

				var orderingAlg = new AmdSymmetricOrdering();
				(int[] permutation, bool oldToNew) = orderingAlg.FindPermutation(order, cscRowIndices, cscColOffsets);
				Assert.True(!oldToNew);
				comparer.AssertEqual(SparsePosDef10by10.MatlabPermutationAMD, permutation);
			});
		}
	}
}
