namespace MGroup.LinearAlgebra.Implementations.NativeWin64.Tests.Reordering
{
	using System;
	using System.Collections;
	using System.Collections.Generic;
	using System.Linq;
	using System.Text;

	using MGroup.LinearAlgebra.Implementations;
	using MGroup.LinearAlgebra.Implementations.Managed;
	using MGroup.LinearAlgebra.Implementations.NativeWin64;
	using MGroup.LinearAlgebra.Matrices;
	using MGroup.LinearAlgebra.Output;
	using MGroup.LinearAlgebra.Reordering;
	using MGroup.LinearAlgebra.Tests;
	using MGroup.LinearAlgebra.Tests.Reordering;
	using MGroup.LinearAlgebra.Tests.TestData;
	using MGroup.LinearAlgebra.Tests.Utilities;
	using Xunit;

	public class SuiteSparseReorderingTests
	{
		private static readonly MatrixComparer comparer = new MatrixComparer(1E-13);
		private static readonly IImplementationProvider provider = new CustomImplementationProvider(
			ManagedBlasProvider.UniqueInstance, ManagedSparseBlasProvider.UniqueInstance,
			ManagedLapackProvider.UniqueInstance, new SuiteSparseReorderingProvider());

		[SkippableFact]
		private static void TestReorderingAmdSuiteSparse()
		{
			Skip.IfNot(TestSettings.TestSuiteSparse, TestSettings.MessageWhenSkippingSuiteSparse);
			AmdSymmetricOrderingTests.TestFindPermutationGivenPattern(provider);
		}

		[SkippableFact]
		private static void TestReorderingCamdSuiteSparse()
		{
			Skip.IfNot(TestSettings.TestSuiteSparse, TestSettings.MessageWhenSkippingSuiteSparse);

			int n = SparsePosDef10by10.Order;
			var pattern = SparsityPatternSymmetric.CreateFromDense(Matrix.CreateFromArray(SparsePosDef10by10.Matrix));
			var orderingAlg = new OrderingCamdSuiteSparse();
			(int[] permutation, ReorderingStatistics stats) =
				orderingAlg.FindPermutation(pattern, SparsePosDef10by10.ConstraintsCAMD);

			var originalDiagonal = new double[n];
			var permutedDiagonal = new double[n];
			for (int i = 0; i < n; ++i)
			{
				originalDiagonal[i] = SparsePosDef10by10.Matrix[i, i];
			}

			for (int i = 0; i < n; ++i)
			{
				permutedDiagonal[i] = originalDiagonal[permutation[i]];
			}

			var writer = new Array1DWriter();
			Console.Write("Permutation (new-to-old): ");
			writer.WriteToConsole(permutation);
			Console.Write("Original diagonal: ");
			writer.WriteToConsole(originalDiagonal);
			Console.Write("Permuted diagonal: ");
			writer.WriteToConsole(permutedDiagonal);

			comparer.AssertEqual(SparsePosDef10by10.PermutationCAMD, permutation);
		}
	}
}
