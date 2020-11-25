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
	/// Eigenvalues and optionally left and right eigenvectors of a real square matrix. The matrix can be nonsymmetric and is
	/// stored in full column major format. Uses Lapack.
	/// </summary>
	public class NonSymmetricEigensystemFull : IEigensystem
	{
		private NonSymmetricEigensystemFull(int order, double[] eigenvaluesReal, double[] eigenvaluesImaginary,
			double[] eigenvectorsLeft, double[] eigenvectorsRight)
		{
			this.Order = order;
			this.EigenvaluesReal = Vector.CreateFromArray(eigenvaluesReal);
			this.EigenvaluesImaginary = Vector.CreateFromArray(eigenvaluesImaginary);

			if (eigenvectorsLeft == null)
			{
				this.EigenvectorsLeft = null;
			}
			else
			{
				this.EigenvectorsLeft = Matrix.CreateFromArray(eigenvectorsLeft, order, order);
			}

			if (eigenvectorsRight == null)
			{
				this.EigenvectorsRight = null;
			}
			else
			{
				this.EigenvectorsRight = Matrix.CreateFromArray(eigenvectorsRight, order, order);
			}
		}

		/// <summary>
		/// See <see cref="IEigensystem.EigenvaluesReal"/>. 
		/// There are <see cref="Order"/> eigenvalues, but their real parts may be repeated.
		/// </summary>
		public Vector EigenvaluesReal { get; }

		/// <summary>
		/// See <see cref="IEigensystem.EigenvaluesImaginary"/>.
		/// There are <see cref="Order"/> eigenvalues. For complex eigenvalues, the conjugate pairs will be stored one after the
		/// other and their correspoding real parts in <see cref="EigenvaluesReal"/> will be identical.
		/// </summary>
		public Vector EigenvaluesImaginary { get; } = null;

		/// <summary>
		/// See <see cref="IEigensystem.EigenvectorsLeft"/>.
		/// </summary>
		public Matrix EigenvectorsLeft { get; }

		/// <summary>
		/// See <see cref="IEigensystem.EigenvectorsRight"/>.
		/// </summary>
		public Matrix EigenvectorsRight { get; }

		/// <summary>
		/// See <see cref="IEigensystem.Order"/>.
		/// </summary>
		public int Order { get; }

		/// <summary>
		/// Calculates the eigenvalues and optionally the left and right eigenvectors of a square matrix. Does not check if the
		/// matrix is square.
		/// </summary>
		/// <param name="order">The number of rows/columns of the original square matrix.</param>
		/// <param name="matrix">
		/// The original matrix in full column major format. Will be overwritten.
		/// </param>
		/// <param name="calcLeftEigenvectors">
		/// If true, left eigenvectors will be computed, in addition to eigenvalues and possibly right eigenvectors.
		/// Else left eigenvectors will not be computed.
		/// </param>
		/// <param name="calcRightEigenvectors">
		/// If true, right eigenvectors will be computed, in addition to eigenvalues and possibly left eigenvectors.
		/// Else right eigenvectors will not be computed.
		/// </param>
		/// <returns>An object holding the eigensystem of the matrix.</returns>
		public static NonSymmetricEigensystemFull Create(
			int order, double[] matrix, bool calcLeftEigenvectors, bool calcRightEigenvectors)
		{
			// Prepare LAPACK input
			int leadingDimA = order;

			// There as many eigenvalues as the order of the matrix, but some the complex eigenvalues come in conjugate pairs.
			var eigenvaluesReal = new double[order];
			var eigenvaluesImaginary = new double[order];

			EigensystemJob jobLeft;
			double[] eigenvectorsLeft;
			int leadingDimEigenvectorsLeft;
			if (calcLeftEigenvectors)
			{
				jobLeft = EigensystemJob.EigenvaluesAndEigenVectors;
				eigenvectorsLeft = new double[order * order];
				leadingDimEigenvectorsLeft = order;
			}
			else
			{
				jobLeft = EigensystemJob.OnlyEigenvalues;
				eigenvectorsLeft = new double[1];
				leadingDimEigenvectorsLeft = 1;
			}

			EigensystemJob jobRight;
			double[] eigenvectorsRight;
			int leadingDimEigenvectorsRight;
			if (calcRightEigenvectors)
			{
				jobRight = EigensystemJob.EigenvaluesAndEigenVectors;
				eigenvectorsRight = new double[order * order];
				leadingDimEigenvectorsRight = order;
			}
			else
			{
				jobRight = EigensystemJob.OnlyEigenvalues;
				eigenvectorsRight = new double[1];
				leadingDimEigenvectorsRight = 1;
			}

			// Call Lapack
			LapackEigensystems.Dgeev(jobLeft, jobRight, order, matrix, 0, leadingDimA, 
				eigenvaluesReal, 0, eigenvaluesImaginary, 0,
				eigenvectorsLeft, 0, leadingDimEigenvectorsLeft,
				eigenvectorsRight, 0, leadingDimEigenvectorsRight);

			// Repack LAPACK output into usable classes
			if (!calcLeftEigenvectors) eigenvectorsLeft = null;
			if (!calcRightEigenvectors) eigenvectorsRight = null;
			return new NonSymmetricEigensystemFull(
				order, eigenvaluesReal, eigenvaluesImaginary, eigenvectorsLeft, eigenvectorsRight);
		}
	}
}
