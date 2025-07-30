using System;
using System.Collections.Generic;
using System.Text;

using MGroup.LinearAlgebra.Exceptions;
using MGroup.LinearAlgebra.Matrices;
using MGroup.LinearAlgebra.Vectors;

namespace MGroup.LinearAlgebra.Distributed.Overlapping
{
    public static class CompatibilityUtilities
    {
		public static DistributedOverlappingVector CastToDistributed(Vectors.IVectorView vector)
		{
			if (vector is DistributedOverlappingVector casted)
			{
				return casted;
			}
			else
			{
				throw new NonMatchingFormatException($"Cannot perform the required operation, because the vector is not in " +
					$"{typeof(DistributedOverlappingVector)} format, but in {vector.GetType()} format.");
			}
		}

		public static DistributedOverlappingMatrix<TMatrix> CastToDistributed<TMatrix>(Matrices.IMatrixView matrix)
			where TMatrix : class, IMatrix
		{
			if (matrix is DistributedOverlappingMatrix<TMatrix> casted)
			{
				return casted;
			}
			else
			{
				throw new NonMatchingFormatException($"Cannot perform the required operation, because the matrix is not in " +
					$"{typeof(DistributedOverlappingMatrix<TMatrix>)} format, but in {matrix.GetType()} format.");
			}
		}

		public static void CheckSameFormat(DistributedOverlappingVector vector1, DistributedOverlappingVector vector2)
		{
			if (!vector2.Indexer.IsCompatibleWith(vector1.Indexer))
			{
				throw new NonMatchingFormatException(
					$"The 2 vectors have different formats, as defined by their indexers ({vector1.Indexer.GetType()})");
			}
		}

		public static void CheckSameFormat<TMatrix>(
			DistributedOverlappingMatrix<TMatrix> matrix1, DistributedOverlappingMatrix<TMatrix> matrix2)
			where TMatrix : class, IMatrix
		{
			if (!matrix1.Indexer.IsCompatibleWith(matrix2.Indexer))
			{
				throw new NonMatchingFormatException(
					$"The 2 matrices have different formats, as defined by their indexers ({matrix1.Indexer.GetType()})");
			}
		}

		public static void CheckSameFormat<TMatrix>(
			DistributedOverlappingMatrix<TMatrix> matrix, DistributedOverlappingVector vector)
			where TMatrix : class, IMatrix
		{
			if (!matrix.Indexer.IsCompatibleWith(vector.Indexer))
			{
				throw new NonMatchingFormatException(
					$"The matrix and the vector have different formats, as defined by their indexers ({matrix.Indexer.GetType()})");
			}
		}

		public static void CheckSameFormat(DistributedOverlappingTransformation matrix, DistributedOverlappingVector vector)
		{
			if (!matrix.Indexer.IsCompatibleWith(vector.Indexer))
			{
				throw new NonMatchingFormatException(
					$"The matrix and the vector have different formats, as defined by their indexers ({matrix.Indexer.GetType()})");
			}
		}
	}
}
