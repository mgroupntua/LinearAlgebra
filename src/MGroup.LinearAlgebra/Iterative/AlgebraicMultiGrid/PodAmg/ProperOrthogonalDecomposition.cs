namespace MGroup.LinearAlgebra.Iterative.AlgebraicMultiGrid.PodAmg
{
	using System;
	using System.Collections.Generic;
	using System.Text;

	using MGroup.LinearAlgebra.Eigensystems;
	using MGroup.LinearAlgebra.Matrices;
	using MGroup.LinearAlgebra.Vectors;

	public class ProperOrthogonalDecomposition
	{
		private readonly bool _keepOnlyNonZeroEigenvalues;
		private readonly double _zeroEigenvalueTolerance;

		public ProperOrthogonalDecomposition(bool keepOnlyNonZeroEigenvalues, double zeroEigenvalueTolerance=1E-10)
		{
			_keepOnlyNonZeroEigenvalues = keepOnlyNonZeroEigenvalues;
			_zeroEigenvalueTolerance = zeroEigenvalueTolerance;
		}

		/// <summary>
		/// Performs POD analysis and returns the principal components of a samples set
		/// </summary>
		/// <param name="numSampleVectors">
		/// The number of sample vectors. Must be equal to the number of columns in <paramref name="sampleVectors"/> and &gt;= 2.
		/// </param>
		/// <param name="sampleVectors">
		/// A (d x n) matrix, which contains n vectors of length d. Each of these n vectors corresponds to one sample, time-step, 
		/// etc. In general n &lt; d, which is used here for increased performance.
		/// </param>
		/// <param name="numPrincipalComponents">
		/// How many principal components to keep. Depending on the configuration of this object, if some eigenvectors 
		/// correspond to zero eigenvalues, they will be discarded and fewer total eigenvectors will be returned.
		/// </param>
		/// <returns></returns>
		public Matrix CalculatePrincipalComponents(int numSampleVectors, Matrix sampleVectors, int numPrincipalComponents) 
		{ 
			if (sampleVectors.NumColumns != numSampleVectors)
			{
				throw new ArgumentException("The matrix containing the sample vectors must have " +
					$"{numSampleVectors} columns but was ({sampleVectors.NumRows}, {sampleVectors.NumColumns}).");
			}
			if (numSampleVectors < 2) 
			{
				throw new ArgumentException("There must be at least 2 vectors (columns) in the input matrix.");
			}
			if (numPrincipalComponents > sampleVectors.NumColumns)
			{
				throw new ArgumentException(
					"Cannot request more principal components than the number of vectors (columns) in the input matrix.");
			}

			// If d > n (usual case): Psi = eigenvectors of U^T*U and Phi = U * Psi
			// If d <= n: Phi = eigenvectors of U*U^T.
			if (sampleVectors.NumRows > sampleVectors.NumColumns)
			{
				Matrix correlation = sampleVectors.MultiplyRight(sampleVectors, transposeThis: true, transposeOther: false);
				var svd = SingularValueDecomposition.Calculate(correlation);

				int numComponentsToKeep = CountPrincipalComponentsToKeep(numPrincipalComponents, svd.SingularValues);
				Matrix principalComponents = sampleVectors * svd.SingularVectors;
				return principalComponents.GetSubmatrix(0, principalComponents.NumRows, 0, numComponentsToKeep); //TODO: discard the unneeded vectors earlier.
			}
			else
			{
				throw new NotImplementedException();
				//Matrix correlation = sampleVectors.MultiplyRight(sampleVectors, transposeThis: false, transposeOther: true);
				//var svd = SingularValueDecomposition.Calculate(correlation);
				//Matrix principalComponents = svd.SingularVectors;
				//return principalComponents.GetSubmatrix(0, principalComponents.NumRows, 0, numPrincipalComponents); //TODO: discard the unneeded vectors earlier.
			}
		}

		private int CountPrincipalComponentsToKeep(int numComponentsRequested, Vector eigenvaluesDescending)
		{
			if (_keepOnlyNonZeroEigenvalues)
			{
				int numComponentsToKeep = 0;
				for (int i = 0; i < numComponentsRequested; ++i)
				{
					if (Math.Abs(eigenvaluesDescending[i]) <= _zeroEigenvalueTolerance) // Only keep eigenvectors of non-zero eigenvalues
					{
						break;
					}
					++numComponentsToKeep;
				}
				return numComponentsToKeep;
			}
			else
			{
				return numComponentsRequested;
			}
		}
	}
}
