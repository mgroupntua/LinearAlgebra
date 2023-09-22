//using MGroup.LinearAlgebra.Iterative;
//using MGroup.LinearAlgebra.Iterative.AlgebraicMultiGrid;
//using MGroup.LinearAlgebra.Iterative.ConjugateGradient;
//using MGroup.LinearAlgebra.Iterative.ConjugateGradient;
//using MGroup.LinearAlgebra.Iterative.Termination;
//using MGroup.LinearAlgebra.Matrices;
//using MGroup.LinearAlgebra.Tests.TestData;
//using MGroup.LinearAlgebra.Tests.Utilities;
//using MGroup.LinearAlgebra.Vectors;
//using Xunit;

//namespace MGroup.LinearAlgebra.Tests.Iterative
//{
//	using LinearAlgebra.Iterative.PreconditionedConjugateGradient;

//	/// <summary>
//    /// Tests for <see cref="AlgebraicMultiGrid"/>.
//    /// Authors: Gerasimos Sotiropoulos 
//    /// </summary>
//    public static class AMGTests
//	{
//        private static readonly MatrixComparer comparer = new MatrixComparer(1E-5);

        
//        [Theory]
//        [MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
//        private static void TestSparseSystem(LinearAlgebraProviderChoice providers)
//        {
//            TestSettings.RunMultiproviderTest(providers, delegate ()
//            {
//                var A = Matrix.CreateFromArray(SparsePosDef10by10.Matrix);
//                var b = Vector.CreateFromArray(SparsePosDef10by10.Rhs);
//                var xExpected = Vector.CreateFromArray(SparsePosDef10by10.Lhs);

//                var builder = new AlgebraicMultiGrid.Builder();
//                builder.ResidualTolerance = 1E-7;
//                builder.MaxIterationsProvider = new PercentageMaxIterationsProvider(1.0);
//                var cg = builder.Build();
//                var xComputed = Vector.CreateZero(A.NumRows);
//                IterativeStatistics stats = cg.Solve(A, b, xComputed);
//                comparer.AssertEqual(xExpected, xComputed);
//            });
//        }
//    }
//}
