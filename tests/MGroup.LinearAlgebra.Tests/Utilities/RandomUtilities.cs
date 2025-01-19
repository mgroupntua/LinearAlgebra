namespace MGroup.LinearAlgebra.Tests.Utilities
{
	using System;
	using System.Collections.Generic;
	using System.Text;

	using MGroup.LinearAlgebra.Matrices.Builders;
	using MGroup.LinearAlgebra.Vectors;

	/// <summary>
	/// Utility methods for generating random matrices and vectors.
	/// </summary>
	public static class RandomUtilities
	{
		public static DokSymmetric CreateRandomMatrix(int order, double nonZeroChance)
		{
			var rand = new Random();
			var dok = DokSymmetric.CreateEmpty(order);
			for (int j = 0; j < order; ++j)
			{
				for (int i = 0; i <= j; ++i)
				{
					if (rand.NextDouble() <= nonZeroChance)
					{
						dok[i, j] = rand.NextDouble();
					}
				}
			}
			return dok;
		}

		public static DokRowMajor CreateRandomSparseMatrix(int numRows, int numCols, double nonZeroChance)
		{
			var rand = new Random();
			var dok = DokRowMajor.CreateEmpty(numRows, numCols);
			for (int i = 0; i < numRows; ++i)
			{
				for (int j = 0; j < numCols; ++j)
				{
					if (rand.NextDouble() <= nonZeroChance)
					{
						dok[i, j] = rand.NextDouble();
					}
				}
			}
			return dok;
		}

		public static Vector CreateRandomVector(int length)
		{
			var rand = new Random();
			var vector = new double[length];
			for (int i = 0; i < length; ++i)
			{
				vector[i] = rand.NextDouble();
			}
			return Vector.CreateFromArray(vector, false);
		}
	}
}
