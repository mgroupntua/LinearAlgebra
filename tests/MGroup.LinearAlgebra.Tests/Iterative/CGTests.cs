using MGroup.LinearAlgebra.Iterative;
using MGroup.LinearAlgebra.Iterative.ConjugateGradient;
using MGroup.LinearAlgebra.Iterative.ConjugateGradient;
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
    public static class CGTests
    {
        private static readonly MatrixComparer comparer = new MatrixComparer(1E-5);

        [Theory]
        [MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
        private static void TestDenseSystem(LinearAlgebraProviderChoice providers)
        {
            TestSettings.RunMultiproviderTest(providers, delegate ()
            {
                var A = Matrix.CreateFromArray(SymmPosDef10by10.Matrix);
                var b = Vector.CreateFromArray(SymmPosDef10by10.Rhs);
                var xExpected = Vector.CreateFromArray(SymmPosDef10by10.Lhs);

                var builder = new CGAlgorithm.Builder();
                builder.ResidualTolerance = 1E-7;
                builder.MaxIterationsProvider = new PercentageMaxIterationsProvider(1.0);
                var cg = builder.Build();
                var xComputed = Vector.CreateZero(A.NumRows);
                IterativeStatistics stats = cg.Solve(A, b, xComputed, true);
                comparer.AssertEqual(xExpected, xComputed);
            });
        }

        [Theory]
        [MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
        private static void TestSparseSystem(LinearAlgebraProviderChoice providers)
        {
            TestSettings.RunMultiproviderTest(providers, delegate ()
            {
                var A = Matrix.CreateFromArray(SparsePosDef10by10.Matrix);
                var b = Vector.CreateFromArray(SparsePosDef10by10.Rhs);
                var xExpected = Vector.CreateFromArray(SparsePosDef10by10.Lhs);

                var builder = new CGAlgorithm.Builder();
                builder.ResidualTolerance = 1E-7;
                builder.MaxIterationsProvider = new PercentageMaxIterationsProvider(1.0);
                var cg = builder.Build();
                var xComputed = Vector.CreateZero(A.NumRows);
                IterativeStatistics stats = cg.Solve(A, b, xComputed, true);
                comparer.AssertEqual(xExpected, xComputed);
            });
        }
    }
}
