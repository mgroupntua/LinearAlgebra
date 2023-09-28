namespace MGroup.LinearAlgebra.Tests.Eigensystems
{
	using System;
	using System.Collections.Generic;
	using System.Text;

	using MGroup.LinearAlgebra.Commons;
	using MGroup.LinearAlgebra.Matrices;
	using MGroup.LinearAlgebra;
	using MGroup.LinearAlgebra.Vectors;
	using MGroup.LinearAlgebra.Tests.Utilities;
	using System.Linq;

	//TODO: Extend this to also support complex eigenvalues, left eigenvectors, and eigenvalues with multiplicity.
	public class SpectralDecomposition
	{
		private readonly double tolerance;
		private readonly bool useRightEigenvectors;
		private readonly Vector eigenvaluesReal;
		private readonly Matrix eigenvectorsRight;

		private readonly SortedSet<double> eigenvalues;
		private readonly Dictionary<double, Vector> rightEigenvectorsOfEigenvalues; // this assumes only one eigenvector per eigenvalue, which is not the case for higher multiplicities.

		public SpectralDecomposition(Vector eigenvaluesReal, Matrix eigenvectorsRight, double tolerance) 
		{
			this.tolerance = tolerance;
			this.useRightEigenvectors = eigenvectorsRight != null;
			this.eigenvaluesReal = eigenvaluesReal;
			this.eigenvectorsRight = eigenvectorsRight;

			var eigenvalueComparer = new EigenvalueComparer(new ValueComparer(tolerance));
			this.eigenvalues = new SortedSet<double>(eigenvalueComparer);
			this.rightEigenvectorsOfEigenvalues = new Dictionary<double, Vector>();

			int n = eigenvaluesReal.Length;
			for (int i = 0; i < n; i++) 
			{
				var lambda = eigenvaluesReal[i]; //TODO: also take care of duplicate eigenvalues
				this.eigenvalues.Add(lambda);
				if (useRightEigenvectors)
				{
					rightEigenvectorsOfEigenvalues[lambda] = eigenvectorsRight.GetColumn(i);
				}
			}
		}

		public bool CanRecomposeOriginalMatrix(IMatrixView originalMatrix)
		{
			if (!useRightEigenvectors)
			{
				return false;
			}

			int n = originalMatrix.NumRows;
			Matrix invX = eigenvectorsRight.Invert();
			//Matrix invX = eigenvectorsRight.Transpose(); // Only for spd matrices.
			Matrix reconstructed = eigenvectorsRight * Matrix.CreateFromDiagonal(n, n, eigenvaluesReal.RawData) * invX;

			var comparer = new MatrixComparer(tolerance);
			return comparer.AreEqual(originalMatrix, reconstructed);
		}

		public bool IsEquivalent(SpectralDecomposition other)
		{
			bool sameEigenvalues = HasSameEigenvalues(other);
			if (!sameEigenvalues)
			{
				return false;
			}

			if (this.useRightEigenvectors != other.useRightEigenvectors)
			{
				return false;
			}
			if (this.rightEigenvectorsOfEigenvalues.Count != other.rightEigenvectorsOfEigenvalues.Count)
			{
				return false;
			}

			Dictionary<double, double> valueMap = MapEigenvaluesThisToOther(other);
			foreach (double thisEigenvalue in this.rightEigenvectorsOfEigenvalues.Keys)
			{
				Vector thisEigenvector = this.rightEigenvectorsOfEigenvalues[thisEigenvalue];
				Vector otherEigenvector = other.rightEigenvectorsOfEigenvalues[valueMap[thisEigenvalue]];
				if (!AreLinearlyDependentVectors(thisEigenvector, otherEigenvector)) 
				{
					return false;
				}
			}
			return true;
		}

		public bool IsIdentical(SpectralDecomposition other)
		{
			var comparer = new MatrixComparer(tolerance);
			if (!comparer.AreEqual(this.eigenvaluesReal, other.eigenvaluesReal))
			{
				return false;
			}
			
			if (this.useRightEigenvectors != other.useRightEigenvectors)
			{
				return false;
			}
			if (this.useRightEigenvectors)
			{
				if (!comparer.AreEqual(this.eigenvectorsRight, other.eigenvectorsRight))
				{
					return false;
				}
			}

			return true;
		}

		private bool AreLinearlyDependentVectors(Vector thisEigenvector, Vector otherEigenvector)
		{
			// Two eigenvectors corresponding to the same eigenvalue (of multiplicity=1) may not be identical, but must be
			// linearly dependent.
			int n = thisEigenvector.Length;
			if (otherEigenvector.Length != n) 
			{
				return false;
			}

			// Find the coefficient between the 2 vectors, examining only the first non-zero entry
			double coeff = 0;
			for (int i = 0; i < n; i++)
			{
				if (Math.Abs(thisEigenvector[i]) > tolerance)
				{
					coeff = otherEigenvector[i] / thisEigenvector[i];
				}
			}
			if (coeff == 0) 
			{
				throw new NotImplementedException("One eigenvector is the zero vector.");
			}


			// Check that corresponding entries are multiples with the same coefficient
			var valueComparer = new ValueComparer(tolerance);
			for (int i = 0; i < n; i++) 
			{ 
				double a = thisEigenvector[i];
				double b = otherEigenvector[i];
				if (Math.Abs(a) <= tolerance)
				{
					if (Math.Abs(b) > tolerance)
					{
						return false;
					}
				}
				else
				{
					if (!valueComparer.AreEqual(b / a, coeff))
					{
						return false;
					}
				}
			}
			return true;
		}

		private bool HasSameEigenvalues(SpectralDecomposition other)
		{
			if (this.eigenvalues.Count != other.eigenvalues.Count)
			{
				return false;
			}

			foreach (double lambda in other.eigenvalues)
			{
				if (!this.eigenvalues.Contains(lambda))
				{
					return false;
				}
			}

			return true;
		}

		/// <summary>
		/// Only call this after making sure that the eigenvalues are indeed (almost) identical
		/// </summary>
		private Dictionary<double, double> MapEigenvaluesThisToOther(SpectralDecomposition other)
		{
			double[] thisEigenValues = this.eigenvalues.ToArray();
			double[] otherEigenValues = other.eigenvalues.ToArray();
			var map = new Dictionary<double, double>();
			for (int i = 0; i < thisEigenValues.Length; i++)
			{
				map[thisEigenValues[i]] = otherEigenValues[i];
			}
			return map;
		}

		private class EigenvalueComparer : IComparer<double>
		{
			private readonly ValueComparer valueComparer;

			public EigenvalueComparer(ValueComparer valueComparer)
			{
				this.valueComparer = valueComparer;
			}

			public int Compare(double x, double y)
			{
				if (valueComparer.AreEqual(x, y))
				{
					return 0;
				}
				else
				{
					return Math.Sign(x - y);
				}
			}
		}
	}
}
