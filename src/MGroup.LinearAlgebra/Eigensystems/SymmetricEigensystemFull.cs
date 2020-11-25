using System;
using System.Collections.Generic;
using System.Text;
using MGroup.LinearAlgebra.Matrices;
using MGroup.LinearAlgebra.Providers;
using MGroup.LinearAlgebra.Vectors;
using static MGroup.LinearAlgebra.LibrarySettings;

namespace MGroup.LinearAlgebra.Eigensystems
{
	/// <summary>
	/// Eigenvalues and optionally eigenvectors of a symmetric matrix stored in full column major format. Only the upper
	/// triangle or lower part of the original matrix is modified. The other part of the original matrix is still stored in the 
	/// array but it is ignored. Uses Lapack.
	/// </summary>
	public class SymmetricEigensystemFull : IEigensystem
	{
		private SymmetricEigensystemFull(int order, Vector eigenvaluesReal, Matrix eigenvectorsRight)
		{
			this.Order = order;
			this.EigenvaluesReal = eigenvaluesReal;
			this.EigenvectorsRight = eigenvectorsRight;
		}

		/// <summary>
		/// See <see cref="IEigensystem.EigenvaluesReal"/>. 
		/// Symmetric matrices have <see cref="Order"/> real eigenvalues, stored here in ascending order.
		/// </summary>
		public Vector EigenvaluesReal { get; }

		/// <summary>
		/// See <see cref="IEigensystem.EigenvaluesImaginary"/>. 
		/// Null since symmetric matrices have real eigenvalues.
		/// </summary>
		public Vector EigenvaluesImaginary { get; } = null;

		/// <summary>
		/// See <see cref="IEigensystem.EigenvectorsRight"/>.
		/// </summary>
		public Matrix EigenvectorsRight { get; }

		/// <summary>
		/// See <see cref="IEigensystem.EigenvectorsLeft"/>.
		/// Null since they are the transpose of <see cref="EigenvectorsRight"/>.
		/// </summary>
		public Matrix EigenvectorsLeft { get; } = null;

		/// <summary>
		/// See <see cref="IEigensystem.Order"/>.
		/// </summary>
		public int Order { get; }

		/// <summary>
		/// Calculates the eigenvalues and optionally the eigenvectors of a symmetric matrix. Does not check if the matrix is
		/// symmetric.
		/// </summary>
		/// <param name="order">The number of rows/columns of the original symmetric matrix.</param>
		/// <param name="matrix">
		/// The original matrix in full column major format. Will be overwritten.
		/// </param>
		/// <param name="calcEigenvectors">
		/// If true, both eigenvalues and eigenvectors will be computed. Else only eigenvalues will be computed.
		/// </param>
		/// <returns>An object holding the eigensystem of the matrix.</returns>
		public static SymmetricEigensystemFull Create(int order, double[] matrix, bool calcEigenvectors)
		{
			EigensystemJob job = calcEigenvectors ? EigensystemJob.EigenvaluesAndEigenVectors : EigensystemJob.OnlyEigenvalues;
			int leadingDimA = order;
			var eigenvalues = new double[order]; // symmetric matrices have as many eigenvalues as their size
			LapackEigensystems.Dsyev(job, StoredTriangle.Upper, order, matrix, 0, leadingDimA, eigenvalues, 0);

			if (calcEigenvectors)
			{
				// The original matrix is overwritten by the eigenvectors
				var eigenvectors = Matrix.CreateFromArray(matrix, order, order);
				return new SymmetricEigensystemFull(order, Vector.CreateFromArray(eigenvalues), eigenvectors);
			}
			else
			{
				return new SymmetricEigensystemFull(order, Vector.CreateFromArray(eigenvalues), null);
			}
		}
	}
}
