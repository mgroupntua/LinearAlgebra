namespace MGroup.LinearAlgebra.Implementations.Managed
{
	using System;

	/// <summary>
	/// Provides managed C# implementations of the linear algebra operations defined by <see cref="IBlasProvider"/>.
	/// Ports FORTRAN BLAS routines, when available. Otherwise, custom C# implementations are used instead.
	/// </summary>
	public partial class ManagedBlasProvider : IBlasProvider
	{
		internal enum Diagonal { Regular, Unit, Zero };

		public static ManagedBlasProvider UniqueInstance { get; } = new ManagedBlasProvider();

		private ManagedBlasProvider() { } // private constructor for singleton pattern

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
