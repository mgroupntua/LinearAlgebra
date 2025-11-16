namespace MGroup.LinearAlgebra.Tests.Utilities
{
	using System;
	using System.Collections.Generic;
	using System.Linq;
	using System.Text;
	using System.Threading.Tasks;

	using MGroup.LinearAlgebra.Matrices;
	using MGroup.LinearAlgebra.Vectors;

	internal class MatrixVectorMultiplicationChecker
	{
		protected MatrixComparer comparer;

		internal MatrixVectorMultiplicationChecker()
		{
			this.comparer = new MatrixComparer(1E-13);
			CreateLhsVectorFunc = x => Vector.CreateFromArray(x, true);
			CreateZeroRhsVectorFunc = n => Vector.CreateZero(n);
			CreateRandomVectorFunc = n =>
			{
				int seed = 13;
				var rng = new Random(seed);
				var result = Vector.CreateZero(n);
				for (int i = 0; i < n; i++)
				{
					result[i] = rng.NextDouble();
				}
				return result;
			};
		}

		internal Func<double[], IReadOnlyVector> CreateLhsVectorFunc { get; set; }

		internal Func<int, IVector> CreateZeroRhsVectorFunc { get; set; }

		internal Func<int, IVector> CreateRandomVectorFunc { get; set; }

		internal double Tolerance
		{
			set => comparer = new MatrixComparer(1E-13);
		}

		internal virtual void CheckMultiplication(IReadOnlyMatrix matrix, double[] lhsVector, double[] rhsVectorExpected, 
			bool transposeMatrix)
		{
			var lhs = CreateLhsVectorFunc(lhsVector);
			IVector rhs = matrix.Multiply(lhs, transposeMatrix);
			comparer.AssertEqual(rhsVectorExpected, rhs);
		}

		internal virtual void CheckMultiplicationIntoResult(IReadOnlyMatrix matrix, double[] lhsVector, double[] rhsVectorExpected, 
			bool transposeMatrix)
		{
			var lhs = CreateLhsVectorFunc(lhsVector);
			IVector rhs = CreateRandomVectorFunc(rhsVectorExpected.Length);
			matrix.MultiplyIntoResult(lhs, rhs, transposeMatrix);
			comparer.AssertEqual(rhsVectorExpected, rhs);
		}

		internal virtual void CheckAllMultiplications(IReadOnlyMatrix matrix, double[] lhsVector, double[] rhsVectorExpected,
			bool transposeMatrix)
		{
			CheckMultiplication(matrix, lhsVector, rhsVectorExpected, transposeMatrix);
			CheckMultiplicationIntoResult(matrix, lhsVector, rhsVectorExpected, transposeMatrix);
		}
	}
}
