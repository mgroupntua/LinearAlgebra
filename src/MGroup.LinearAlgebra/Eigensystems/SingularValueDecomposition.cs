namespace MGroup.LinearAlgebra.Eigensystems
{
	using System;
	using System.Collections;
	using System.Collections.Generic;
	using System.Linq;
	using System.Text;

	using MGroup.LinearAlgebra.Commons;
	using MGroup.LinearAlgebra.Matrices;
	using MGroup.LinearAlgebra.Vectors;

	/// <summary>
	/// See the following links for LAPACK implementations:
	/// https://www.intel.com/content/www/us/en/docs/onemkl/code-samples-lapack/2023-1/gesvd-function.html,
	/// https://www.intel.com/content/www/us/en/docs/onemkl/code-samples-lapack/2023-1/gesdd-function.html,
	/// https://www.intel.com/content/www/us/en/docs/onemkl/code-samples-lapack/2023-1/lapacke-dgesvd-example-c-column.html
	/// https://www.ibm.com/docs/en/essl/6.2?topic=llss-sgesvd-dgesvd-cgesvd-zgesvd-sgesdd-dgesdd-cgesdd-zgesdd-singular-value-decomposition-general-matrix
	/// </summary>
	public class SingularValueDecomposition //TODO: The singular vectors are not identical to numpy's
	{
		private SingularValueDecomposition(Vector singularValues, Matrix singularVectors) 
		{
			SingularValues = singularValues;
			SingularVectors = singularVectors;
		}

		/// <summary>
		/// Non negative, in descending order
		/// </summary>
		public Vector SingularValues { get; }

		public Matrix SingularVectors { get; }

		//TODO: Should all factorization classes have as input the containers such as Matrix, instead of raw arrays?
		//		It would make them easier to consume. E.g. POD should not have to concern itself with the raw arrays of a general
		//		or packed matrix. It should only see Matrix or SymmetricMatrix.
		public static SingularValueDecomposition Calculate(Matrix matrix) 
		{
			Preconditions.CheckSquare(matrix);
			int n = matrix.NumRows;
			var singularValues = new double[n];
			var singularVectors = new double[n, n];
			DenseStrategies.SVD(matrix, singularValues, singularVectors);
			//return new SingularValueDecomposition(Vector.CreateFromArray(singularValues), Matrix.CreateFromArray(singularVectors));

			//TODO: sorting and storage formats should be handled inside this SVD for more efficiency
			(Vector valuesSorted, Matrix vectorsSorted) = SpectralUtilities.SortSingularValues(
				Vector.CreateFromArray(singularValues), Matrix.CreateFromArray(singularVectors), descending: true);
			if (valuesSorted[n - 1] < 0)
			{
				//TODO: add tolerance, dedicated exception class and perhaps make this test optional.
				throw new Exception("Negative singular value found");
			}

			return new SingularValueDecomposition(valuesSorted, vectorsSorted);
		}
	}
}
