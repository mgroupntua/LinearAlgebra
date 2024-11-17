namespace MGroup.LinearAlgebra.Tests.SchurComplements
{
	using System;
	using System.Collections.Generic;
	using System.Linq;
	using System.Text;
	using System.Threading.Tasks;

	using MGroup.LinearAlgebra.Matrices;
	using MGroup.LinearAlgebra.SchurComplements;
	using MGroup.LinearAlgebra.Tests.TestData;
	using MGroup.LinearAlgebra.Tests.Utilities;
	using MGroup.LinearAlgebra.Triangulation;

	using Xunit;

	public static class SchurComplementPckCsrSymCscTests
	{
		[Theory]
		[InlineData("A", 0, 1E-14)]
		[InlineData("A", 1, 1E-14)]
		[InlineData("B", 0, 1E-14)]
		[InlineData("B", 1, 1E-14)]
		public static void TestCalcSchurComplement(string exampleLetter, int indexGroup, double comparisonTolerance)
		{
			SchurComplementExample example;
			if (exampleLetter == "A")
			{
				example = SchurComplementExample.CreateExampleSymmetricA();
			}
			else if (exampleLetter == "B")
			{
				example = SchurComplementExample.CreateExampleSymmetricB();
			}
			else
			{
				throw new NotImplementedException();
			}

			SymmetricMatrix A00;
			CsrMatrix A01;
			SymmetricCscMatrix A11;
			double[,] expectedS11;

			int n = example.MatrixOrder;
			if (indexGroup == 0)
			{
				int n0 = example.Indices1.Length;
				int n1 = example.Indices0.Length;
				A00 = SymmetricMatrix.CreateFromPackedColumnMajorArray(example.Submatrix11PackedUpper, copyArray: true);
				A01 = CsrMatrix.CreateFromArrays(n0, n1, example.Submatrix10Csr.values, example.Submatrix10Csr.colIndices,
					example.Submatrix10Csr.rowOffsets, true);
				A11 = SymmetricCscMatrix.CreateFromArrays(n1, example.Submatrix00Csc.values, example.Submatrix00Csc.rowIndices,
					example.Submatrix00Csc.colOffsets, true);
				expectedS11 = example.SchurOfSubmatrix00;
			}
			else if (indexGroup == 1)
			{
				int n0 = example.Indices0.Length;
				int n1 = example.Indices1.Length;
				A00 = SymmetricMatrix.CreateFromPackedColumnMajorArray(example.Submatrix00PackedUpper, copyArray: true);
				A01 = CsrMatrix.CreateFromArrays(n0, n1, example.Submatrix01Csr.values, example.Submatrix01Csr.colIndices,
					example.Submatrix01Csr.rowOffsets, true);
				A11 = SymmetricCscMatrix.CreateFromArrays(n1, example.Submatrix11Csc.values, example.Submatrix11Csc.rowIndices,
					example.Submatrix11Csc.colOffsets, true);
				expectedS11 = example.SchurOfSubmatrix11;
			}
			else
			{
				throw new ArgumentOutOfRangeException("Valid values for index group are 0, 1", nameof(indexGroup));
			}

			CalcAndTestSchurComplement(A00, A01, A11, comparisonTolerance, expectedS11);
		}

		private static void CalcAndTestSchurComplement(SymmetricMatrix submatrix00, CsrMatrix submatrix01,
			SymmetricCscMatrix submatrix11, double comparisonTolerance, double[,] expectedSchur11)
		{
			var comparer = new MatrixComparer(comparisonTolerance);
			using var inverseA11 = new CholeskyCSparseNet();
			inverseA11.Factorize(submatrix11);

			// Test the method that returns a new instance for the Schur complement
			SymmetricMatrix computedS11 = SchurComplementPckCsrSymCsc.CalcSchurComplement(submatrix00, submatrix01, inverseA11);
			comparer.AssertEqual(expectedSchur11, computedS11);

			// Test the method that overwrites an existing matrix with the Schur complement
			computedS11 = SymmetricMatrix.CreateZero(submatrix00.NumRows);
			SchurComplementPckCsrSymCsc.CalcSchurComplement(submatrix00, submatrix01, inverseA11, computedS11);
			comparer.AssertEqual(expectedSchur11, computedS11);
		}
	}
}
