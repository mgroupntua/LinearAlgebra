using MGroup.LinearAlgebra.Triangulation;
using MGroup.LinearAlgebra.Matrices;
using MGroup.LinearAlgebra.Matrices.Builders;
using MGroup.LinearAlgebra.Tests.TestData;
using MGroup.LinearAlgebra.Tests.Utilities;
using MGroup.LinearAlgebra.Vectors;
using Xunit;
using MGroup.LinearAlgebra.Implementations.Managed;
using MGroup.LinearAlgebra.Implementations;
using MGroup.LinearAlgebra.Implementations.NativeWin64.Triangulation;

namespace MGroup.LinearAlgebra.Tests.Triangulation
{
	/// <summary>
	/// Tests for <see cref="CholeskyCSparseNet"/>.
	/// </summary>
	public static class CholeskySymmetricCscTests
	{
		private static readonly MatrixComparer comparer = new MatrixComparer(1E-13);

		[Theory]
		[MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
		private static void TestSystemSolution1(IImplementationProvider provider)
		{
			TestSettings.RunMultiproviderTest(provider, delegate ()
			{
				int order = SparsePosDef10by10.Order;
				var skyline = SkylineMatrix.CreateFromArrays(order, SparsePosDef10by10.SkylineValues,
					SparsePosDef10by10.SkylineDiagOffsets, true, true);
				var dok = DokSymmetric.CreateFromSparseMatrix(skyline);
				var b = Vector.CreateFromArray(SparsePosDef10by10.Rhs);
				var xExpected = Vector.CreateFromArray(SparsePosDef10by10.Lhs);

				(double[] cscValues, int[] cscRowIndices, int[] cscColOffsets) = dok.BuildSymmetricCscArrays(true);

				using ICholeskySymmetricCsc factorization = provider.CreateSymmetricCscTriangulation(superNodal: true);
				factorization.Factorize(order, cscValues.Length, cscValues, cscRowIndices, cscColOffsets);
				Vector xComputed = factorization.SolveLinearSystem(b);
				comparer.AssertEqual(xExpected, xComputed);
			});
		}

		[Theory]
		[MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
		private static void TestSystemSolution2(IImplementationProvider provider)
		{
			TestSettings.RunMultiproviderTest(provider, delegate ()
			{
				// Define linear system
				var rhs = Vector.CreateFromArray(new double[] { 6.0, 14.0, 11.0, 12.0 });
				var solutionExpected = Vector.CreateFromArray(new double[] { 1.0, 1.0, 1.0, 1.0 });
				var matrixDOK = DokSymmetric.CreateEmpty(4);
				matrixDOK[0, 0] = 4.0; matrixDOK[0, 2] = 2.0;
				matrixDOK[1, 1] = 10.0; matrixDOK[1, 2] = 1.0; matrixDOK[1, 3] = 3.0;
				matrixDOK[2, 2] = 8.0;
				matrixDOK[3, 3] = 9.0;
				SymmetricCscMatrix matrixCSC = matrixDOK.BuildSymmetricCscMatrix(true);

				//const int n = 4;
				//const int nnz = 7;
				//int[] colOffsets = new int[n + 1] { 0, 1, 2, 5, nnz };
				//int[] rowIndices = new int[nnz] { 0, 1, 0, 1, 2, 1, 3 };
				//double[] values = new double[nnz] { 4.0, 10.0, 2.0, 1.0, 8.0, 3.0, 9.0 };
				//SymmetricCSC matrixCSC = new SymmetricCSC(values, rowIndices, colOffsets, false);

				using ICholeskySymmetricCsc factorization = provider.CreateSymmetricCscTriangulation(superNodal: true);
				factorization.Factorize(matrixCSC);
				Vector solution = factorization.SolveLinearSystem(rhs);
				comparer.AssertEqual(solutionExpected, solution);
			});
		}
	}
}
