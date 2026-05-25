namespace MGroup.LinearAlgebra.Implementations.Managed.Custom
{
	using System;

	/// <summary>
	/// Provides managed C# implementations of the linear algebra operations defined by <see cref="IBlasProvider"/>.
	/// </summary>
	public partial class CustomBlasProvider : IBlasProvider
	{
		internal enum Diagonal { Regular, Unit, Zero };

		public static CustomBlasProvider UniqueInstance { get; } = new CustomBlasProvider();

		private CustomBlasProvider() { } // private constructor for singleton pattern

		private static bool UseUpperImplementation(StoredTriangle uplo, TransposeMatrix transA)
		{
			//if (transA == TransposeMatrix.ConjugateTranspose)
			//    throw new ArgumentException("Cannot use conjugate transpose operations for double matrices and vectors.");

			if (uplo == StoredTriangle.Upper && transA == TransposeMatrix.NoTranspose) return true;
			if (uplo == StoredTriangle.Upper && transA == TransposeMatrix.Transpose) return false;
			if (uplo == StoredTriangle.Lower && transA == TransposeMatrix.NoTranspose) return true;
			if (uplo == StoredTriangle.Lower && transA == TransposeMatrix.Transpose) return false;
			throw new Exception("This code should not have been reached");
		}
	}
}
