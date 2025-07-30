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
		}

		internal Func<double[], IVectorView> CreateLhsVectorFunc { get; set; }

		internal Func<int, IVector> CreateZeroRhsVectorFunc { get; set; }

		internal double Tolerance
		{
			set => comparer = new MatrixComparer(1E-13);
		}

		internal virtual void CheckMultiplication(IMatrixView matrix, double[] lhsVector, double[] rhsVectorExpected, 
			bool transposeMatrix)
		{
			var lhs = CreateLhsVectorFunc(lhsVector);
			IVector rhs = matrix.Multiply(lhs, transposeMatrix);
			comparer.AssertEqual(rhsVectorExpected, rhs);
		}

		internal virtual void CheckMultiplicationIntoResult(IMatrixView matrix, double[] lhsVector, double[] rhsVectorExpected, 
			bool transposeMatrix)
		{
			var lhs = CreateLhsVectorFunc(lhsVector);
			IVector rhs = CreateZeroRhsVectorFunc(rhsVectorExpected.Length);
			matrix.MultiplyIntoResult(lhs, rhs, transposeMatrix);
			comparer.AssertEqual(rhsVectorExpected, rhs);
		}

		internal virtual void CheckAllMultiplications(IMatrixView matrix, double[] lhsVector, double[] rhsVectorExpected,
			bool transposeMatrix)
		{
			CheckMultiplication(matrix, lhsVector, rhsVectorExpected, transposeMatrix);
			CheckMultiplicationIntoResult(matrix, lhsVector, rhsVectorExpected, transposeMatrix);
		}
	}
}
