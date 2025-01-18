//TODO: Implement IIndexable2D to allow easy output.
//TODO: Allow other orderings, as I do in CholeskySuiteSparse
//TODO: CSparse.NET also provides matrix update and downdate operations.
//TODO: Improve error checking
//TODO: expose internal factorized arrays
namespace MGroup.LinearAlgebra.Triangulation
{
	using System;

	using CSparse;
	using CSparse.Double;
	using CSparse.Double.Factorization;

	using MGroup.LinearAlgebra.Exceptions;
	using MGroup.LinearAlgebra.Matrices;
	using MGroup.LinearAlgebra.Vectors;

	/// <summary>
	/// LU factorization of a sparse square matrix using the CSparse.NET library. The original matrix must be in
	/// Compressed Sparse Columns format.
	/// </summary>
	public class LUCSparseNet : ILUCscFactorization
	{
		public const double DefaultPivotTolelance = 0.1; // between 0.0, 1.0. TODO: Find a good default

		private SparseLU factorization;

		public int NumColumns => Order;

		/// <summary>
		/// The number of non-zero entries (and explicitly stored zeros) in the explicitly stored lower and upper triangular
		/// factors after LU factorization.
		/// </summary>
		public int NumNonZerosUpper => factorization.NonZerosCount;

		public int NumRows => Order;

		/// <summary>
		/// The number of rows/columns of the square matrix.
		/// </summary>
		public int Order { get; private set; }

		/// <summary>
		/// The internal data containing the factorization.
		/// </summary>
		public SparseLU RawData => factorization;

		public void Dispose() { } // Do nothing. This is purely managed code.

		/// <inheritdoc/>
		public void Factorize(
			int order, int numNonZeros, double[] cscValues, int[] cscRowIndices, int[] cscColOffsets, double pivotTolerance)
		{
			try
			{
				var matrixCSparse = new SparseMatrix(order, order, cscValues, cscRowIndices, cscColOffsets);
				var factorization = SparseLU.Create(matrixCSparse, ColumnOrdering.Natural, pivotTolerance);

				this.Order = order;
				this.factorization = factorization;
			}
			catch (Exception ex) //TODO: how can I make sure this exception was thrown because of an indefinite matrix?
			{
				throw new SingularMatrixException(ex.Message);
			}
		}

		/// <inheritdoc/>
		public void Factorize(CscMatrix matrix, double pivotTolerance)
			=> Factorize(matrix.NumColumns, matrix.NumNonZeros, matrix.RawValues, matrix.RawRowIndices, matrix.RawColOffsets,
				pivotTolerance);

		/// <summary>
		/// See <see cref="ITriangulation.CalcDeterminant"/>.
		/// </summary>
		public double CalcDeterminant()
		{
			throw new NotImplementedException();
		}

		/// <summary>
		/// See <see cref="ITriangulation.SolveLinearSystem(Vector, Vector)"/>.
		/// </summary>
		public void SolveLinearSystem(Vectors.Vector rhs, Vectors.Vector solution)
			=> factorization.Solve(rhs.RawData, solution.RawData);
	}
}
