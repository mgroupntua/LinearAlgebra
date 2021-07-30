using System;
using System.Collections.Generic;
using System.Text;

namespace MGroup.LinearAlgebra.Tests.TestData
{
    internal static class RandomMatrices
    {
        private const int defaultSeed = 273;

        internal static double[,] CreateRandomMatrix(int numRows, int numCols, int seed = defaultSeed)
        {
            var rng = new Random(seed);
            var matrix = new double[numRows, numCols];
            for (int i = 0; i < numRows; ++i)
            {
                for (int j = 0; j < numCols; ++j) matrix[i, j] = rng.NextDouble();
            }
            return matrix;
        }

        internal static double[,] CreateRandomMatrix(int numRows, int numCols, double min, double max, int seed = defaultSeed)
        {
            var rng = new Random(seed);
            var matrix = new double[numRows, numCols];
            for (int i = 0; i < numRows; ++i)
            {
                for (int j = 0; j < numCols; ++j) matrix[i, j] = min + (max - min) * rng.NextDouble();
            }
            return matrix;
        }

        internal static double[] CreateRandomVector(int length, int seed = defaultSeed)
        {
            var rng = new Random(seed);
            var vector = new double[length];
            for (int i = 0; i < length; ++i) vector[i] = rng.NextDouble();
            return vector;
        }

        internal static double[] CreateRandomVector(int length, double min, double max, int seed = defaultSeed)
        {
            var rng = new Random(seed);
            var vector = new double[length];
            for (int i = 0; i < length; ++i) vector[i] = min + (max - min) * rng.NextDouble();
            return vector;
        }
    }
}
