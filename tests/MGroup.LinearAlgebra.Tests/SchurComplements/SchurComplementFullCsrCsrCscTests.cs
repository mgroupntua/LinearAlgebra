namespace MGroup.LinearAlgebra.Tests.SchurComplements
{
	using System;
	using System.Collections.Generic;
	using System.Linq;
	using System.Text;
	using System.Threading.Tasks;

	using MGroup.LinearAlgebra.Implementations;
	using MGroup.LinearAlgebra.Implementations.Managed;
	using MGroup.LinearAlgebra.Matrices;
	using MGroup.LinearAlgebra.SchurComplements;
	using MGroup.LinearAlgebra.Tests.TestData;
	using MGroup.LinearAlgebra.Tests.Utilities;
	using MGroup.LinearAlgebra.Triangulation;

	using Xunit;

	public static class SchurComplementFullCsrCsrCscTests
	{
		[Theory]
		[InlineData("A", false, 0, 1E-14)]
		[InlineData("A", false, 1, 1E-14)]
		[InlineData("A", true, 0, 1E-14)]
		[InlineData("A", true, 1, 1E-14)]
		[InlineData("B", false, 0, 1E-14)]
		[InlineData("B", false, 1, 1E-14)]
		[InlineData("B", true, 0, 1E-14)]
		[InlineData("B", true, 1, 1E-14)]
		public static void TestCalcSchurComplement(string exampleLetter, bool symmetric, int indexGroup, 
			double comparisonTolerance)
		{
			SchurComplementExample example;
			if (exampleLetter == "A")
			{
				if (symmetric)
				{
					example = SchurComplementExample.CreateExampleSymmetricA();
				}
				else
				{
					example = SchurComplementExample.CreateExampleNonSymmetricA();
				}
			}
			else if (exampleLetter == "B")
			{
				if (symmetric)
				{
					example = SchurComplementExample.CreateExampleSymmetricB();
				}
				else
				{
					example = SchurComplementExample.CreateExampleNonSymmetricB();
				}
			}
			else
			{
				throw new NotImplementedException();
			}

			Matrix A00;
			CsrMatrix A01, A10;
			CscMatrix A11;
			double[,] expectedS11;

			int n = example.MatrixOrder;
			if (indexGroup == 0)
			{
				int n0 = example.Indices1.Length;
				int n1 = example.Indices0.Length;
				A00 = Matrix.CreateFromArray(example.Submatrix11FullColMajor, n0, n0, true);
				A01 = CsrMatrix.CreateFromArrays(n0, n1, example.Submatrix10Csr.values, example.Submatrix10Csr.colIndices,
					example.Submatrix10Csr.rowOffsets, true);
				A10 = CsrMatrix.CreateFromArrays(n1, n0, example.Submatrix01Csr.values, example.Submatrix01Csr.colIndices,
					example.Submatrix01Csr.rowOffsets, true);
				A11 = CscMatrix.CreateFromArrays(n1, n1, example.Submatrix00Csc.values, example.Submatrix00Csc.rowIndices,
					example.Submatrix00Csc.colOffsets, true);
				expectedS11 = example.SchurOfSubmatrix00;
			}
			else if (indexGroup == 1)
			{
				int n0 = example.Indices0.Length;
				int n1 = example.Indices1.Length;
				A00 = Matrix.CreateFromArray(example.Submatrix00FullColMajor, n0, n0, true);
				A01 = CsrMatrix.CreateFromArrays(n0, n1, example.Submatrix01Csr.values, example.Submatrix01Csr.colIndices,
					example.Submatrix01Csr.rowOffsets, true);
				A10 = CsrMatrix.CreateFromArrays(n1, n0, example.Submatrix10Csr.values, example.Submatrix10Csr.colIndices,
					example.Submatrix10Csr.rowOffsets, true);
				A11 = CscMatrix.CreateFromArrays(n1, n1, example.Submatrix11Csc.values, example.Submatrix11Csc.rowIndices,
					example.Submatrix11Csc.colOffsets, true);
				expectedS11 = example.SchurOfSubmatrix11;
			}
			else
			{
				throw new ArgumentOutOfRangeException("Valid values for index group are 0, 1", nameof(indexGroup));
			}

			CalcAndTestSchurComplement(A00, A01, A10, A11, comparisonTolerance, expectedS11);
		}

		private static void CalcAndTestSchurComplement(Matrix submatrix00, CsrMatrix submatrix01, CsrMatrix submatrix10,
			CscMatrix submatrix11, double comparisonTolerance, double[,] expectedSchur11)
		{
			var comparer = new MatrixComparer(comparisonTolerance);

			// Factorize
			IImplementationProvider provider = new ManagedSequentialImplementationProvider();
			ILUCscFactorization inverseA11 = provider.CreateLUTriangulation();
			inverseA11.Factorize(submatrix11, 1E-7);

			// Test the method that returns a new instance for the Schur complement
			Matrix computedS11 = SchurComplementFullCsrCsrCsc.CalcSchurComplement(submatrix00, submatrix01, submatrix10, inverseA11);
			comparer.AssertEqual(expectedSchur11, computedS11);

			// Test the method that overwrites an existing matrix with the Schur complement
			computedS11 = Matrix.CreateZero(submatrix00.NumRows, submatrix00.NumColumns);
			SchurComplementFullCsrCsrCsc.CalcSchurComplement(submatrix00, submatrix01, submatrix10, inverseA11, computedS11);
			comparer.AssertEqual(expectedSchur11, computedS11);
		}
	}
}
