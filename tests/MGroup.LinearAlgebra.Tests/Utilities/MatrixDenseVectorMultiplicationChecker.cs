namespace MGroup.LinearAlgebra.Tests.Utilities
{
	using System;
	using System.Collections.Generic;
	using System.Linq;
	using System.Text;
	using System.Threading.Tasks;

	using MGroup.LinearAlgebra.Matrices;
	using MGroup.LinearAlgebra.Vectors;

	internal class MatrixDenseVectorMultiplicationChecker : MatrixVectorMultiplicationChecker
	{
		internal delegate Vector MultiplyDenseVector(IReadOnlyMatrix matrix, Vector lhs, bool transposeMatrix);

		internal delegate void MultiplyDenseVectorIntoResult(IReadOnlyMatrix matrix, Vector lhs, Vector rhs, bool transposeMatrix);

		private readonly MultiplyDenseVector multiplyVectorFunc;
		private readonly MultiplyDenseVectorIntoResult multiplyVectorIntoResultFunc;

		internal MatrixDenseVectorMultiplicationChecker(MultiplyDenseVector multiplyVectorFunc,
			MultiplyDenseVectorIntoResult multiplyVectorIntoResultFunc)
			: base()
		{
			this.multiplyVectorFunc = multiplyVectorFunc;
			this.multiplyVectorIntoResultFunc = multiplyVectorIntoResultFunc;
		}

		internal void CheckMultiplicationDense(IReadOnlyMatrix matrix, double[] lhsVector, double[] rhsVectorExpected, 
			bool transposeMatrix)
		{
			var lhs = Vector.CreateFromArray(lhsVector, true);
			Vector rhs = multiplyVectorFunc(matrix, lhs, transposeMatrix);
			comparer.AssertEqual(rhsVectorExpected, rhs);
		}

		internal void CheckMultiplicationIntoResultDense(IReadOnlyMatrix matrix, double[] lhsVector, double[] rhsVectorExpected, 
			bool transposeMatrix)
		{
			var lhs = Vector.CreateFromArray(lhsVector, true);
			var rhs = Vector.CreateWithValue(rhsVectorExpected.Length, -3.33);
			multiplyVectorIntoResultFunc(matrix, lhs, rhs, transposeMatrix);
			comparer.AssertEqual(rhsVectorExpected, rhs);
		}

		internal override void CheckAllMultiplications(IReadOnlyMatrix matrix, double[] lhsVector, double[] rhsVectorExpected,
			bool transposeMatrix)
		{
			CheckMultiplication(matrix, lhsVector, rhsVectorExpected, transposeMatrix);
			CheckMultiplicationIntoResult(matrix, lhsVector, rhsVectorExpected, transposeMatrix);
			CheckMultiplicationDense(matrix, lhsVector, rhsVectorExpected, transposeMatrix);
			CheckMultiplicationIntoResultDense(matrix, lhsVector, rhsVectorExpected, transposeMatrix);
		}
	}
}
