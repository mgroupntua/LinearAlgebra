namespace MGroup.LinearAlgebra.Matrices
{
	using System;
	using System.Collections.Generic;
	using System.Text;
	using MGroup.LinearAlgebra.Commons;
	using MGroup.LinearAlgebra.Exceptions;
	using MGroup.LinearAlgebra.Vectors;

	public interface ΙDiagonalAccessible : IIndexable2D
	{
		/// <summary>
		/// Returns a <see cref="Vector"/> with the entries of the matrix's main diagonal.
		/// </summary>
		/// <exception cref="NonMatchingDimensionsException">Thrown if the matrix is not square.</exception>
		public Vector GetDiagonal() => Vector.CreateFromArray(this.GetDiagonalAsArray(), false);

		/// <summary>
		/// Returns an array with the entries of the matrix's main diagonal.
		/// </summary>
		/// <exception cref="NonMatchingDimensionsException">Thrown if the matrix is not square.</exception>
		public double[] GetDiagonalAsArray()
		{
			Preconditions.CheckSquare(this);
			double[] diag = new double[this.NumRows];
			for (int i = 0; i < this.NumRows; ++i)
			{
				diag[i] = this[i, i];
			}
			return diag;
		}
	}
}
