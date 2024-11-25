namespace MGroup.LinearAlgebra.Tests.Reordering
{
	using System;
	using System.IO;
	using System.Reflection;

	using MGroup.LinearAlgebra.Matrices;
	using MGroup.LinearAlgebra.Output;
	using MGroup.LinearAlgebra.Reordering;
	using MGroup.LinearAlgebra.Tests.TestData;
	using MGroup.LinearAlgebra.Tests.Utilities;

	using Xunit;

	/// <summary>
	/// Tests for <see cref="SparsityPatternSymmetric"/>.
	/// </summary>
	public static class SparsityPatternSymmetricTests
	{
		[Fact]
		private static void TestConnectIndices()
		{
			int n = GlobalMatrixAssembly.GlobalOrder;
			var dense = Matrix.CreateFromArray(GlobalMatrixAssembly.GlobalMatrix);
			var pattern = SparsityPatternSymmetric.CreateEmpty(n);
			pattern.ConnectIndices(GlobalMatrixAssembly.GlobalIndices1, true);
			pattern.ConnectIndices(GlobalMatrixAssembly.GlobalIndices2, true);
			pattern.ConnectIndices(GlobalMatrixAssembly.GlobalIndices3, true);

			for (int i = 0; i < n; ++i)
			{
				for (int j = 0; j < n; ++j)
				{
					bool denseHasZero = dense[i, j] == 0.0;
					bool patternHasZero = !pattern.IsNonZero(i, j);
					Assert.True(patternHasZero == denseHasZero);
				}
			}
		}
	}
}
