using MGroup.LinearAlgebra.Commons;
using MGroup.LinearAlgebra.Exceptions;
using MGroup.LinearAlgebra.Matrices;
using MGroup.LinearAlgebra.Tests.TestData;
using MGroup.LinearAlgebra.Tests.Utilities;
using MGroup.LinearAlgebra.Vectors;
using Xunit;

namespace MGroup.LinearAlgebra.Tests.Matrices
{
    /// <summary>
    /// Tests for <see cref="TriangularUpper"/>.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class TriangularUpperTests
    {
        private static readonly MatrixComparer comparer = new MatrixComparer(1E-13);

        [Fact]
        private static void TestArrayCopy()
        {
            // invertible
            var A1 = TriangularUpper.CreateFromArray(UpperInvertible10by10.Matrix);
            comparer.AssertEqual(UpperInvertible10by10.Matrix, A1.CopyToArray2D());

            // singular
            var A2 = TriangularUpper.CreateFromArray(UpperSingular10by10.Matrix);
            comparer.AssertEqual(UpperSingular10by10.Matrix, A2.CopyToArray2D());
        }

        [Fact]
        private static void TestClear()
        {
            var zero = Matrix.CreateZero(UpperSingular10by10.Order, UpperSingular10by10.Order);
            var matrix = TriangularUpper.CreateFromArray(UpperSingular10by10.Matrix);
            matrix.Clear();
            comparer.AssertEqual(zero, matrix);
        }

		[SkippableFact(typeof(PerformanceBottleneckException))]
        private static void TestGetColumn()
		{
            var matrix = TriangularUpper.CreateFromArray(UpperInvertible10by10.Matrix);
            for (int j = 0; j < UpperInvertible10by10.Order; ++j)
            {
                Vector colExpected = DenseStrategies.GetColumn(matrix, j);
                Vector colComputed = matrix.GetColumn(j);
                comparer.AssertEqual(colExpected, colComputed);
            }
        }

		[SkippableFact(typeof(PerformanceBottleneckException))]
		private static void TestGetRow()
        {
            var matrix = TriangularUpper.CreateFromArray(UpperInvertible10by10.Matrix);
            for (int i = 0; i < UpperInvertible10by10.Order; ++i)
            {
                Vector rowExpected = DenseStrategies.GetRow(matrix, i);
                Vector rowComputed = matrix.GetRow(i);
                comparer.AssertEqual(rowExpected, rowComputed);
            }
        }

        [Theory]
        [MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
        private static void TestMatrixVectorMultiplication(LinearAlgebraProviderChoice providers)
        {
            TestSettings.RunMultiproviderTest(providers, delegate ()
            {
                // invertible
                var A1 = TriangularUpper.CreateFromArray(UpperInvertible10by10.Matrix);
                var x1 = Vector.CreateFromArray(UpperInvertible10by10.Lhs);
                var b1Expected = Vector.CreateFromArray(UpperInvertible10by10.Rhs);
                Vector b1Computed = A1.Multiply(x1);
                comparer.AssertEqual(b1Expected, b1Computed);

                // singular
                var A2 = TriangularUpper.CreateFromArray(UpperSingular10by10.Matrix);
                var x2 = Vector.CreateFromArray(UpperSingular10by10.Lhs);
                var b2Expected = Vector.CreateFromArray(UpperSingular10by10.Rhs);
                Vector b2Computed = A2.Multiply(x1);
                comparer.AssertEqual(b2Expected, b2Computed);
            });
        }

        [Theory]
        [MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
        private static void TestMatrixVectorMultiplicationIntoResult(LinearAlgebraProviderChoice providers)
        {
            // The result vectors will first be set to some non zero values to make sure that the result overwrites 
            // them instead of being added to them.

            TestSettings.RunMultiproviderTest(providers, delegate ()
            {
                var A = TriangularUpper.CreateFromArray(UpperInvertible10by10.Matrix);
                var x = Vector.CreateFromArray(UpperInvertible10by10.Lhs);
                var bExpected = Vector.CreateFromArray(UpperInvertible10by10.Rhs);
                Vector bComputed = Vector.CreateWithValue(A.NumRows, 1.0);
                A.MultiplyIntoResult(x, bComputed, false);
                comparer.AssertEqual(bExpected, bComputed);
            });
        }

        [Theory]
        [MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
        private static void TestSystemSolution(LinearAlgebraProviderChoice providers)
        {
            TestSettings.RunMultiproviderTest(providers, delegate ()
            {
                // invertible
                var A1 = TriangularUpper.CreateFromArray(UpperInvertible10by10.Matrix);
                var b1 = Vector.CreateFromArray(UpperInvertible10by10.Rhs);
                var x1Expected = Vector.CreateFromArray(UpperInvertible10by10.Lhs);
                Vector x1Computed = A1.SolveLinearSystem(b1);
                comparer.AssertEqual(x1Expected, x1Computed);

                // singular
                var A2 = TriangularUpper.CreateFromArray(UpperSingular10by10.Matrix);
                var b2 = Vector.CreateFromArray(UpperSingular10by10.Rhs);
                var x2Expected = Vector.CreateFromArray(UpperSingular10by10.Lhs);
                Vector x2Computed = A2.SolveLinearSystem(b2);
                Assert.False(comparer.AreEqual(x2Expected, x2Computed));

                // invertible - solve transposed (forward substitution)
                Matrix A3 = Matrix.CreateFromArray(UpperInvertible10by10.Matrix).Invert().Transpose();
                Vector x3Expected = A3 * b1;
                Vector x3Computed = A1.SolveLinearSystem(b1, true);
                comparer.AssertEqual(x3Expected, x3Computed);
            });
        }

        [Fact]
        private static void TestTransposition()
        {
            // invertible
            var A1 = TriangularUpper.CreateFromArray(UpperInvertible10by10.Matrix);
            var A1TransposeExpected = MatrixOperations.Transpose(UpperInvertible10by10.Matrix);
            TriangularLower A1TransposeComputed = A1.Transpose(false);
            comparer.AssertEqual(A1TransposeExpected, A1TransposeComputed.CopyToArray2D());

            // singular
            var A2 = TriangularUpper.CreateFromArray(UpperSingular10by10.Matrix);
            var A2TransposeExpected = MatrixOperations.Transpose(UpperSingular10by10.Matrix);
            TriangularLower A2TransposeComputed = A2.Transpose(false);
            comparer.AssertEqual(A2TransposeExpected, A2TransposeComputed.CopyToArray2D());
        }
    }
}
