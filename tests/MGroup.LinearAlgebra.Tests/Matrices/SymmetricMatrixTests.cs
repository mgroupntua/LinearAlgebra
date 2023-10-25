using System.IO;
using System.Runtime.Serialization.Formatters.Binary;
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
	/// Tests for <see cref="SymmetricMatrix"/>.
	/// Authors: Serafeim Bakalakos
	/// </summary>
	public static class SymmetricMatrixTests
	{
		private static readonly MatrixComparer comparer = new MatrixComparer(1E-13);

		[Fact]
		private static void TestArrayCopy()
		{
			// positive definite
			var A1 = SymmetricMatrix.CreateFromArray(SymmPosDef10by10.Matrix);
			comparer.AssertEqual(SymmPosDef10by10.Matrix, A1.CopyToArray2D());

			// singular
			var A2 = SymmetricMatrix.CreateFromArray(SymmSingular10by10.Matrix);
			comparer.AssertEqual(SymmSingular10by10.Matrix, A2.CopyToArray2D());
		}

		[Fact]
		private static void TestClear()
		{
			var zero = Matrix.CreateZero(SymmPosDef10by10.Order, SymmPosDef10by10.Order);
			var matrix = SymmetricMatrix.CreateFromArray(SymmPosDef10by10.Matrix);
			matrix.Clear();
			comparer.AssertEqual(zero, matrix);
		}

		[SkippableFact(typeof(PerformanceBottleneckException))]
		private static void TestGetColumn()
		{
			var matrix = SymmetricMatrix.CreateFromArray(SymmPosDef10by10.Matrix);
			for (int j = 0; j < SymmPosDef10by10.Order; ++j)
			{
				Vector colExpected = DenseStrategies.GetColumn(matrix, j);
				Vector colComputed = matrix.GetColumn(j);
				comparer.AssertEqual(colExpected, colComputed);
			}
		}

		[SkippableFact(typeof(PerformanceBottleneckException))]
		private static void TestGetRow()
		{
			var matrix = SymmetricMatrix.CreateFromArray(SymmPosDef10by10.Matrix);
			for (int i = 0; i < SymmPosDef10by10.Order; ++i)
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
				var A1 = SymmetricMatrix.CreateFromArray(SymmPosDef10by10.Matrix);
				var x1 = Vector.CreateFromArray(SymmPosDef10by10.Lhs);
				var b1Expected = Vector.CreateFromArray(SymmPosDef10by10.Rhs);
				Vector b1Computed = A1.Multiply(x1);
				comparer.AssertEqual(b1Expected, b1Computed);

				// singular
				var A2 = SymmetricMatrix.CreateFromArray(SymmSingular10by10.Matrix);
				var x2 = Vector.CreateFromArray(SymmSingular10by10.Lhs);
				var b2Expected = Vector.CreateFromArray(SymmSingular10by10.Rhs);
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
				var A = SymmetricMatrix.CreateFromArray(SymmPosDef10by10.Matrix);
				var x = Vector.CreateFromArray(SymmPosDef10by10.Lhs);
				var bExpected = Vector.CreateFromArray(SymmPosDef10by10.Rhs);
				Vector bComputed = Vector.CreateWithValue(A.NumRows, 1.0);
				A.MultiplyIntoResult(x, bComputed, false);
				comparer.AssertEqual(bExpected, bComputed);
			});
		}

		[SkippableFact(typeof(PerformanceBottleneckException))]
		private static void TestSerialization()
		{
			var originalMatrix = SymmetricMatrix.CreateFromArray(SymmPosDef10by10.Matrix);
			var formatter = new BinaryFormatter();
			using (var stream = new MemoryStream())
			{
				formatter.Serialize(stream, originalMatrix);
				stream.Seek(0, SeekOrigin.Begin);
				var deserializedMatrix = (SymmetricMatrix)formatter.Deserialize(stream);

				Assert.True(originalMatrix.Equals(deserializedMatrix));
			}
		}

		[Fact]
		private static void TestTransposition()
		{
			// invertible
			var A1 = SymmetricMatrix.CreateFromArray(SymmPosDef10by10.Matrix);
			var A1TransposeExpected = MatrixOperations.Transpose(SymmPosDef10by10.Matrix);
			SymmetricMatrix A1TransposeComputed = A1.Transpose(false);
			comparer.AssertEqual(A1TransposeExpected, A1TransposeComputed.CopyToArray2D());

			// singular
			var A2 = SymmetricMatrix.CreateFromArray(SymmSingular10by10.Matrix);
			var A2TransposeExpected = MatrixOperations.Transpose(SymmSingular10by10.Matrix);
			SymmetricMatrix A2TransposeComputed = A2.Transpose(false);
			comparer.AssertEqual(A2TransposeExpected, A2TransposeComputed.CopyToArray2D());
		}
	}
}
