//TODO: each matrix-vector multiplication will internally cast the vectors, in order to perform efficiently. Could I avoid all
//      that casting by using generics
namespace MGroup.LinearAlgebra.Iterative
{
	using MGroup.LinearAlgebra.Matrices;
	using MGroup.LinearAlgebra.Vectors;

	/// <summary>
	/// Wrapper for a matrix class so that it can be used by iterative algorithms, which operate on
	/// <see cref="ILinearTransformation"/>.
	/// </summary>
	public class ExplicitMatrixTransformation : ILinearTransformation
	{
		private readonly IReadOnlyMatrix matrix;

		/// <summary>
		/// Initializes a new instance of <see cref="ExplicitMatrixTransformation"/> that wraps the provided <paramref name="matrix"/>.
		/// </summary>
		/// <param name="matrix">The matrix that will be multiplied with vectors during the iterative algorithms.</param>
		public ExplicitMatrixTransformation(IReadOnlyMatrix matrix) => this.matrix = matrix;

		/// <summary>
		/// <inheritdoc/>
		/// </summary>
		public int NumColumns => matrix.NumColumns;

		/// <summary>
		/// <inheritdoc/>
		/// </summary>
		public int NumRows => matrix.NumRows;

		/// <summary>
		/// <inheritdoc/>
		/// </summary>
		public void Multiply(IReadOnlyVector lhsVector, IVector rhsVector) => matrix.MultiplyIntoResult(lhsVector, rhsVector, false);
	}
}
