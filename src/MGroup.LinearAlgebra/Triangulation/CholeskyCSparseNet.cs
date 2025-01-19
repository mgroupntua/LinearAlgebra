using System;

using CSparse;
using CSparse.Double;
using CSparse.Double.Factorization;

using MGroup.LinearAlgebra.Commons;
using MGroup.LinearAlgebra.Exceptions;
using MGroup.LinearAlgebra.Matrices;
using MGroup.LinearAlgebra.Vectors;

//TODO: Implement IIndexable2D to allow easy output.
//TODO: Allow other orderings, as I do in CholeskySuiteSparse
//TODO: CSparse.NET also provides matrix update and downdate operations.
//TODO: Improve error checking
namespace MGroup.LinearAlgebra.Triangulation
{
	/// <summary>
	/// Cholesky factorization of a sparse symmetric positive definite matrix using the CSparse.NET library. The original matrix
	/// must be in Compressed Sparse Columns format, with only the upper triangle stored. This class may serve as a managed
	/// alternative to <see cref="CholeskySuiteSparse"/>.
	/// </summary>
	public class CholeskyCSparseNet : ICholeskySymmetricCsc
	{
		private SparseCholesky factorization;

		public int NumColumns => Order;

		/// <summary>
		/// The number of non-zero entries (and explicitly stored zeros) in the explicitly stored upper triangular factor
		/// after Cholesky factorization.
		/// </summary>
		public int NumNonZerosUpper => factorization.NonZerosCount;

		public int NumRows => Order;

		/// <summary>
		/// The number of rows/columns of the square matrix.
		/// </summary>
		public int Order { get; private set; } = -1;

		/// <summary>
		/// The internal data containing the factorization.
		/// </summary>
		public SparseCholesky RawData => factorization;

		/// <summary>
		/// See <see cref="ITriangulation.CalcDeterminant"/>.
		/// </summary>
		public double CalcDeterminant()
		{
			throw new NotImplementedException();
		}

		public void Dispose() { } // Do nothing. This is purely managed code.

		/// <inheritdoc/>
		public void Factorize(int order, int numNonZerosUpper, double[] cscValues, int[] cscRowIndices, int[] cscColOffsets)
		{
			if (factorization != null)
			{
				throw new InvalidOperationException("A factorization already exists.");
			}

			try
			{
				var matrixCSparse = new SparseMatrix(order, order, cscValues, cscRowIndices, cscColOffsets);
				var factorization = SparseCholesky.Create(matrixCSparse, ColumnOrdering.Natural);

				this.Order = order;
				this.factorization = factorization;
			}
			catch (Exception ex) //TODO: how can I make sure this exception was thrown because of an indefinite matrix?
			{
				throw new IndefiniteMatrixException(ex.Message);
			}
		}

		/// <inheritdoc/>
		public void Factorize(SymmetricCscMatrix matrix) => Factorize(
			matrix.NumColumns, matrix.NumNonZerosUpper, matrix.RawValues, matrix.RawRowIndices, matrix.RawColOffsets);

		/// <summary>
		/// See <see cref="ITriangulation.SolveLinearSystem(Vector, Vector)"/>.
		/// </summary>
		public void SolveLinearSystem(Vectors.Vector rhs, Vectors.Vector solution)
		{
			if (factorization == null)
			{
				throw new InvalidOperationException("No matrix has been factorized yet");
			}

			factorization.Solve(rhs.RawData, solution.RawData);
		}

		/// <inheritdoc/>
		public Matrix SolveLinearSystems(Matrix rhsVectors)
		{
			Preconditions.CheckSystemSolutionDimensions(Order, rhsVectors.NumRows);
			var result = Matrix.CreateZero(Order, rhsVectors.NumColumns);
			var x = Vectors.Vector.CreateZero(Order);
			for (int j = 0; j < rhsVectors.NumColumns; j++)
			{
				Vectors.Vector b = rhsVectors.GetColumn(j);
				factorization.Solve(b.RawData, x.RawData);
				result.SetSubcolumn(j, x, 0);
			}

			return result;
		}
	}
}
