namespace MGroup.LinearAlgebra.Tests.SchurComplements.SubmatrixExtractors
{
	using System;

	using MGroup.LinearAlgebra.Matrices;
	using MGroup.LinearAlgebra.SchurComplements.SubmatrixExtractors;
	using MGroup.LinearAlgebra.Tests.TestData;
	using MGroup.LinearAlgebra.Tests.Utilities;

	using Xunit;

	public static class SubmatrixExtractorPckCsrCscSymTests
	{
		[Fact]
		public static void TestClear()
		{
			var extractor = new SubmatrixExtractorPckCsrCscSym();
			var example = SchurComplementExample.CreateExampleSymmetricA();

			int n = example.MatrixOrder;
			(double[] values, int[] rowIndices, int[] colOffsets) A = example.MatrixCscSymmetric;
			var symmCscA = SymmetricCscMatrix.CreateFromArrays(n, A.values, A.rowIndices, A.colOffsets, true);

			int[] indices0 = example.Indices0;
			int[] indices1 = example.Indices1;
			extractor.ExtractSubmatrices(symmCscA, indices0, indices1);
			TestExtractorForExample(extractor, example);

			extractor.Clear();

			Assert.Null(extractor.Submatrix00);
			Assert.Null(extractor.Submatrix01);
			Assert.Null(extractor.Submatrix11);
		}

		[Fact]
		public static void TestExtractSubmatricesMultiple()
		{
			var extractor = new SubmatrixExtractorPckCsrCscSym();
			var example = SchurComplementExample.CreateExampleSymmetricA();

			int n = example.MatrixOrder;
			(double[] values, int[] rowIndices, int[] colOffsets) A = example.MatrixCscSymmetric;
			var symmCsc = SymmetricCscMatrix.CreateFromArrays(n, A.values, A.rowIndices, A.colOffsets, true);

			int[] indices0 = example.Indices0;
			int[] indices1 = example.Indices1;
			extractor.ExtractSubmatrices(symmCsc, indices0, indices1);

			TestExtractorForExample(extractor, example);

			// 2nd matrix 3*A. Pattern is the same
			double factor = 3.0;
			symmCsc.ScaleIntoThis(factor);
			example.ScalingFactor = factor;
			extractor.ExtractSubmatrices(symmCsc, indices0, indices1);
			TestExtractorForExample(extractor, example);

			// 3rd matrix B is completely different
			example = SchurComplementExample.CreateExampleSymmetricB();
			(double[] values, int[] rowIndices, int[] colOffsets) B = example.MatrixCscSymmetric;
			symmCsc = SymmetricCscMatrix.CreateFromArrays(n, B.values, B.rowIndices, B.colOffsets, true);

			// Without first clearing, it will fail (should it though?)
			indices0 = example.Indices0;
			indices1 = example.Indices1;
			Assert.Throws<InvalidOperationException>(() => extractor.ExtractSubmatrices(symmCsc, indices0, indices1));

			extractor.Clear();
			extractor.ExtractSubmatrices(symmCsc, indices0, indices1);

			TestExtractorForExample(extractor, example);
		}

		[Fact]
		public static void TestExtractSubmatricesOnce()
		{
			var extractor = new SubmatrixExtractorPckCsrCscSym();
			var example = SchurComplementExample.CreateExampleSymmetricA();

			int n = example.MatrixOrder;
			(double[] values, int[] rowIndices, int[] colOffsets) A = example.MatrixCscSymmetric;
			var symmCscA = SymmetricCscMatrix.CreateFromArrays(n, A.values, A.rowIndices, A.colOffsets, true);

			int[] indices0 = example.Indices0;
			int[] indices1 = example.Indices1;
			extractor.ExtractSubmatrices(symmCscA, indices0, indices1);

			TestExtractorForExample(extractor, example);
		}

		private static void TestExtractorForExample(SubmatrixExtractorPckCsrCscSym extractor, SchurComplementExample example)
		{
			double[] expected00 = example.Submatrix00PackedUpper;
			(double[] values, int[] colIndices, int[] rowOffsets) expected01 = example.Submatrix01Csr;
			(double[] values, int[] rowIndices, int[] colOffsets) expected11 = example.Submatrix11CscSymmetric;

			var comparer = new MatrixComparer(1E-15);
			comparer.AssertEqual(expected00, extractor.Submatrix00.RawValues);

			comparer.AssertEqual(expected01.values, extractor.Submatrix01.RawValues);
			comparer.AssertEqual(expected01.colIndices, extractor.Submatrix01.RawColIndices);
			comparer.AssertEqual(expected01.rowOffsets, extractor.Submatrix01.RawRowOffsets);

			comparer.AssertEqual(expected11.values, extractor.Submatrix11.RawValues);
			comparer.AssertEqual(expected11.rowIndices, extractor.Submatrix11.RawRowIndices);
			comparer.AssertEqual(expected11.colOffsets, extractor.Submatrix11.RawColOffsets);
		}
	}
}
