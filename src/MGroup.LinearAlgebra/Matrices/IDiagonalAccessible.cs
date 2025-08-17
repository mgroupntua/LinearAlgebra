namespace MGroup.LinearAlgebra.Matrices
{
	using System;
	using System.Collections.Generic;
	using System.Text;
	using MGroup.LinearAlgebra.Exceptions;
	using MGroup.LinearAlgebra.Vectors;

	public interface IDiagonalAccessible : IIndexable2D
	{
		/// <summary>
		/// Returns a <see cref="Vector"/> with the entries of the matrix's main diagonal. The matrix must be square.
		/// </summary>
		/// <exception cref="NonMatchingDimensionsException">Thrown if the matrix is not square.</exception>
		public Vector GetDiagonal();

		/// <summary>
		/// Returns an array with the entries of the matrix's main diagonal.
		/// </summary>
		/// <exception cref="NonMatchingDimensionsException">Thrown if the matrix is not square.</exception>
		public double[] GetDiagonalAsArray();
	}
}
