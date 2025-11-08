using System;
using System.IO;
using System.Linq;
using System.Runtime.Serialization.Formatters.Binary;

using MGroup.LinearAlgebra.Commons;
using MGroup.LinearAlgebra.Implementations;
using MGroup.LinearAlgebra.Tests.Mocking;
using MGroup.LinearAlgebra.Tests.TestData;
using MGroup.LinearAlgebra.Tests.Utilities;
using MGroup.LinearAlgebra.Vectors;
using Xunit;

namespace MGroup.LinearAlgebra.Tests.Vectors
{
	/// <summary>
	/// Tests for <see cref="DefaultVector"/>.
	/// </summary>
	public static class DefaultVector
	{
		private static readonly MatrixComparer comparer = new MatrixComparer(1E-10);

		[Fact]
		internal static void TestAdd()
		{
			var v1 = new MockDefaultVector(TestVectors.Vector1);
			var v2 = new MockDefaultVector(TestVectors.Vector2);
			var expected = Vector.CreateFromArray(TestVectors.Sum);

			// Add()
			comparer.AssertEqual(expected, v1.Add(v2));

			// AddIntoThis()
			var temp = Vector.CreateFromVector(v1);
			temp.AddIntoThis(v2);
			comparer.AssertEqual(expected, temp);
		}

		[Fact]
		internal static void TestAddIntoThisNonContiguously()
		{
			var v1 = new MockDefaultVector([0.0, 100.0, 200.0, 300.0, 400.0, 500.0, 600.0]);
			var v2 = Vector.CreateFromArray([20.0, 30.0, 50.0], true);
			var expected1 = Vector.CreateFromArray([0.0, 100.0, 220.0, 330.0, 400.0, 550.0, 600.0], true);
			v1.AddIntoThisNonContiguouslyFrom([2, 3, 5], v2);
			comparer.AssertEqual(expected1, v1);

			var v3 = new MockDefaultVector([0.0, 100.0, 200.0, 300.0, 400.0, 500.0, 600.0]);
			var v4 = Vector.CreateFromArray([-10.0, 0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0], true);
			var expected2 = Vector.CreateFromArray([0.0, 110.0, 200.0, 300.0, 440.0, 550.0, 660.0], true);
			v3.AddIntoThisNonContiguouslyFrom([1, 4, 5, 6], v4, [2, 5, 6, 7]);
			comparer.AssertEqual(expected2, v3);
		}

		[Fact]
		internal static void TestAddToIndex()
		{
			var v1 = new MockDefaultVector([0.0, 100.0, 200.0, 300.0, 400.0, 500.0, 600.0]);
			v1.AddToIndex(3, 33.33);
			var expected1 = Vector.CreateFromArray([0.0, 100.0, 200.0, 333.33, 400.0, 500.0, 600.0], true);
			comparer.AssertEqual(expected1, v1);
		}

		[Fact]
		internal static void TestAxpy()
		{
			var v1 = new MockDefaultVector(TestVectors.Vector1);
			var v2 = new MockDefaultVector(TestVectors.Vector2);
			var expected = Vector.CreateFromArray(TestVectors.Vector1PlusVector2Times3);

			// Axpy()
			comparer.AssertEqual(expected, v1.Axpy(v2, TestVectors.Scalar2));

			// AxpyIntoThis
			var computed = new MockDefaultVector(TestVectors.Vector1);
			computed.AxpyIntoThis(v2, TestVectors.Scalar2);
			comparer.AssertEqual(expected, computed);
		}

		[Fact]
		internal static void TestAxpySubvector()
		{
			var v1 = new MockDefaultVector([0.0, 100.0, 200.0, 300.0, 400.0, 500.0, 600.0]);
			var v2 = Vector.CreateFromArray([-10.0, 0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0], true);
			var expected = Vector.CreateFromArray([0.0, 100.0, 200.0, 303.0, 404.0, 505.0, 606.0], true);
			v1.AxpySubvectorIntoThis(3, v2, 0.1, 4, 4);
			comparer.AssertEqual(expected, v1);
		}

		[Fact]
		internal static void TestCopy()
		{
			var original = new MockDefaultVector(TestVectors.Vector1);
			var expected = Vector.CreateFromArray(TestVectors.Vector1);
			IVector clone = original.Copy();
			comparer.AssertEqual(expected, clone);
		}

		[Fact]
		internal static void TestCopyFrom()
		{
			var v1 = new MockDefaultVector(TestVectors.Vector1);
			var v2 = new MockDefaultVector(TestVectors.Vector2);
			var expected = Vector.CreateFromArray(TestVectors.Vector2);
			v1.CopyFrom(v2);
			comparer.AssertEqual(expected, v1);
		}

		[Fact]
		internal static void TestCopyNonContiguouslyFrom()
		{
			var v1 = new MockDefaultVector([20.0, 30.0, 50.0]);
			var v2 = Vector.CreateFromArray([0.0, 100.0, 200.0, 300.0, 400.0, 500.0, 600.0], true);
			var expected1 = Vector.CreateFromArray([200.0, 300.0, 500.0], true);
			v1.CopyNonContiguouslyFrom(v2, [2, 3, 5]);
			comparer.AssertEqual(expected1, v1);

			var v3 = new MockDefaultVector([0.0, 100.0, 200.0, 300.0, 400.0, 500.0, 600.0]);
			var v4 = Vector.CreateFromArray([-10.0, 0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0], true);
			var expected2 = Vector.CreateFromArray([0.0, 10.0, 200.0, 300.0, 40.0, 50.0, 60.0], true);
			v3.CopyNonContiguouslyFrom([1, 4, 5, 6], v4, [2, 5, 6, 7]);
			comparer.AssertEqual(expected2, v3);
		}

		[Fact]
		internal static void TestCopySubvector()
		{
			var v1 = new MockDefaultVector([0.0, 100.0, 200.0, 300.0, 400.0, 500.0, 600.0]);
			var v2 = Vector.CreateFromArray([-10.0, 0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0], true);
			var expected = Vector.CreateFromArray([0.0, 100.0, 200.0, 30.0, 40.0, 50.0, 60.0], true);
			v1.CopySubvectorFrom(3, v2, 4, 4);
			comparer.AssertEqual(expected, v1);
		}

		[Fact]
		internal static void TestCopyToArray()
		{
			var original = new MockDefaultVector(TestVectors.Vector1);
			var expected = TestVectors.Vector1;
			double[] clone = original.CopyToArray();
			comparer.AssertEqual(expected, clone);
		}

		[Fact]
		internal static void TestDoEntrywise()
		{
			var v1 = new MockDefaultVector(TestVectors.Vector1);
			var v2 = new MockDefaultVector(TestVectors.Vector2);
			var expected = Vector.CreateFromArray(
				MatrixOperations.LinearCombination(2.5, TestVectors.Vector1, -3.5, TestVectors.Vector2));
			var comparer = new MatrixComparer();

			// DoEntrywise()
			comparer.AssertEqual(expected, v1.DoEntrywise(v2, (x, y) => 2.5 * x - 3.5 * y));

			// DoEntrywiseIntoThis()
			var computed = new MockDefaultVector(TestVectors.Vector1);
			computed.DoEntrywiseIntoThis(v2, (x, y) => 2.5 * x - 3.5 * y);
			comparer.AssertEqual(expected, computed);
		}

		[Fact]
		internal static void TestDoToAllEntries()
		{
			var vector = new MockDefaultVector(TestVectors.Vector1);
			var expected = Vector.CreateFromArray(TestVectors.Vector1Times2);

			// DoToAllEntries()
			comparer.AssertEqual(expected, vector.DoToAllEntries(x => 2.0 * x));

			// DoToAllEntriesIntoThis()
			var computed = new MockDefaultVector(TestVectors.Vector1);
			computed.DoToAllEntriesIntoThis(x => 2.0 * x);
			comparer.AssertEqual(expected, computed);
		}

		[Fact]
		internal static void TestDotProduct()
		{
			var v1 = new MockDefaultVector(TestVectors.Vector1);
			var v2 = new MockDefaultVector(TestVectors.Vector2);

			// DotProduct()
			comparer.AssertEqual(TestVectors.DotProduct, v1.DotProduct(v2));
		}

		[Fact]
		internal static void TestEquals()
		{
			var v1 = new MockDefaultVector(TestVectors.Vector1);
			var v2 = new MockDefaultVector(TestVectors.Vector1);
			var v3 = Vector.CreateFromArray(TestVectors.Vector1);
			var v4 = Vector.CreateFromArray(TestVectors.Vector2);

			Assert.True(v1.Equals(v2));
			Assert.True(v2.Equals(v1));
			Assert.True(v1.Equals(v3));
			Assert.True(v3.Equals(v1));
			Assert.False(v1.Equals(v4));
			Assert.False(v4.Equals(v1));
		}

		[Fact]
		internal static void TestHadamardProduct()
		{
			var v1 = new MockDefaultVector(TestVectors.Vector1);
			var v2 = new MockDefaultVector(TestVectors.Vector2);
			var expected = Vector.CreateFromArray(TestVectors.HadamardProduct);

			// MultiplyPointwise()
			comparer.AssertEqual(expected, v1.MultiplyEntrywise(v2));

			// MultiplyPointwiseIntoThis()
			var computed = new MockDefaultVector(TestVectors.Vector1);
			computed.MultiplyEntrywiseIntoThis(v2);
			comparer.AssertEqual(expected, computed);
		}

		[Fact]
		internal static void TestLinearCombination()
		{
			var v1 = new MockDefaultVector(TestVectors.Vector1);
			var v2 = new MockDefaultVector(TestVectors.Vector2);
			var expected = Vector.CreateFromArray(
				MatrixOperations.LinearCombination(2.5, TestVectors.Vector1, -3.5, TestVectors.Vector2));
			var comparer = new MatrixComparer();

			// LinearCombination()
			comparer.AssertEqual(expected, v1.LinearCombination(2.5, v2, -3.5));

			// LinearCombinationIntoThis()
			var computed = new MockDefaultVector(TestVectors.Vector1);
			computed.LinearCombinationIntoThis(2.5, v2, -3.5);
			comparer.AssertEqual(expected, computed);
		}

		[Fact]
		internal static void TestNorm2()
		{
			var vector = new MockDefaultVector(TestVectors.Vector1);
			comparer.AssertEqual(TestVectors.Norm2OfVector1, vector.Norm2());
		}

		[Fact]
		internal static void TestReduce()
		{
			var vector = new MockDefaultVector(TestVectors.Difference);
			double maxAbsExpected = MatrixOperations.ReduceMaxAbs(TestVectors.Difference);

			double maxAbsComputed = vector.Reduce(0.0,
				(x, amax) => Math.Abs(x) > amax ? Math.Abs(x) : amax,
				(nz, amax) => amax,
				amax => amax);

			Assert.Equal(maxAbsExpected, maxAbsComputed);
		}

		[Fact]
		internal static void TestScale()
		{
			var vector = new MockDefaultVector(TestVectors.Vector1);
			var expected = Vector.CreateFromArray(TestVectors.Vector1Times2);

			// Scale()
			comparer.AssertEqual(expected, vector.Scale(2.0));

			// ScaleIntoThis()
			var computed = new MockDefaultVector(TestVectors.Vector1);
			computed.ScaleIntoThis(2.0);
			comparer.AssertEqual(expected, computed);
		}

		[Fact]
		internal static void TestSet()
		{
			var v1 = new MockDefaultVector([0.0, 100.0, 200.0, 300.0, 400.0, 500.0, 600.0]);
			v1.Set(3, 33.33);
			var expected1 = Vector.CreateFromArray([0.0, 100.0, 200.0, 33.33, 400.0, 500.0, 600.0], true);
			comparer.AssertEqual(expected1, v1);
		}

		[Fact]
		internal static void TestSetAll()
		{
			var v1 = new MockDefaultVector([0.0, 100.0, 200.0, 300.0, 400.0, 500.0, 600.0]);
			v1.SetAll(-1.1);
			var expected1 = Vector.CreateFromArray([-1.1, -1.1, -1.1, -1.1, -1.1, -1.1, -1.1], true);
			comparer.AssertEqual(expected1, v1);
		}

		[Fact]
		internal static void TestSubtract()
		{
			var v1 = new MockDefaultVector(TestVectors.Vector1);
			var v2 = new MockDefaultVector(TestVectors.Vector2);
			var expected = Vector.CreateFromArray(TestVectors.Difference);

			// Subtract()
			comparer.AssertEqual(expected, v1.Subtract(v2));

			// SubtractIntoThis()
			var computed = new MockDefaultVector(TestVectors.Vector1);
			computed.SubtractIntoThis(v2);
			comparer.AssertEqual(expected, computed);
		}
	}
}
