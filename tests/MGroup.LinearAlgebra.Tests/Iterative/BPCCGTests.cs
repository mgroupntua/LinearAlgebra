using MGroup.LinearAlgebra.Iterative;
using MGroup.LinearAlgebra.Iterative.BlockPreconditionedConjugateGradient;
using MGroup.LinearAlgebra.Iterative.Termination;
using MGroup.LinearAlgebra.Matrices;
using MGroup.LinearAlgebra.Tests.TestData;
using MGroup.LinearAlgebra.Tests.Utilities;
using MGroup.LinearAlgebra.Vectors;

using Xunit;

namespace MGroup.LinearAlgebra.Tests.Iterative
{
	/// <summary>
	/// Tests for <see cref="PcgAlgorithm"/>.
	/// Authors: Serafeim Bakalakos
	/// </summary>
	public static class BPCCGTests
    {
        private static readonly MatrixComparer comparer = new MatrixComparer(1E-5);

        [Theory]
        [MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
        private static void TestDenseSystemStatic(LinearAlgebraProviderChoice providers)
        {
            TestSettings.RunMultiproviderTest(providers, delegate ()
            {
                var A = Matrix.CreateFromArray(SymmPosDef10by10.Matrix);
                var b = Vector.CreateFromArray(SymmPosDef10by10.Rhs);
                var xExpected = Vector.CreateFromArray(SymmPosDef10by10.Lhs);
				
				BPCCGAlgorithm.IterativeStatistics stats = BPCCGAlgorithm.Solve(A, null, b, null,
					5, new PercentageMaxIterationsProvider(1), new RegularResidualConvergence(1e-7));
				var xComputed = stats.solution;
				comparer.AssertEqual(xExpected, xComputed);
			});
        }

		[Theory]
		[MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
		private static void TestDenseSystemObject(LinearAlgebraProviderChoice providers)
		{
			TestSettings.RunMultiproviderTest(providers, delegate ()
			{
				var A = Matrix.CreateFromArray(SymmPosDef10by10.Matrix);
				var b = Vector.CreateFromArray(SymmPosDef10by10.Rhs);
				var xExpected = Vector.CreateFromArray(SymmPosDef10by10.Lhs);
				
				// block size of 6, hits instability
				var bpccg = new BPCCGAlgorithm(new RegularResidualConvergence(1e-7), 5);
				bpccg.Solve(A, b);
				var xComputed = bpccg.Solution;
				comparer.AssertEqual(xExpected, xComputed);
			});
		}
		
		[Theory]
        [MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
        private static void TestSparseSystemStatic(LinearAlgebraProviderChoice providers)
        {
            TestSettings.RunMultiproviderTest(providers, delegate ()
            {
                var A = Matrix.CreateFromArray(SparsePosDef10by10.Matrix);
                var b = Vector.CreateFromArray(SparsePosDef10by10.Rhs);
                var xExpected = Vector.CreateFromArray(SparsePosDef10by10.Lhs);

				BPCCGAlgorithm.IterativeStatistics stats = BPCCGAlgorithm.Solve(A, null, b, null,
					6, new PercentageMaxIterationsProvider(1), new RegularResidualConvergence());
				var xComputed = stats.solution;
				comparer.AssertEqual(xExpected, xComputed);				
            });
        }

		[Theory]
		[MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
		private static void TestSparseSystemObject(LinearAlgebraProviderChoice providers)
		{
			TestSettings.RunMultiproviderTest(providers, delegate ()
			{
				var A = Matrix.CreateFromArray(SparsePosDef10by10.Matrix);
				var b = Vector.CreateFromArray(SparsePosDef10by10.Rhs);
				var xExpected = Vector.CreateFromArray(SparsePosDef10by10.Lhs);

				var bpccg = new BPCCGAlgorithm();
				bpccg.Solve(A, b);
				var xComputed = bpccg.Solution;
				comparer.AssertEqual(xExpected, xComputed);
			});
		}
	}
}
