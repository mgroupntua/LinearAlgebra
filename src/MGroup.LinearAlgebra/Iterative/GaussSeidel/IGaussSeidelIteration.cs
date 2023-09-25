namespace MGroup.LinearAlgebra.Iterative.GaussSeidel
{
	using System;
	using System.Collections.Generic;
	using System.Text;
	using MGroup.LinearAlgebra.Vectors;

	public interface IGaussSeidelIteration : IDisposable
	{
		int SystemSize { get; }

		void Initialize();

		void GaussSeidelBackwardIteration(IVectorView rhsVector, IVector lhsVector);

		void GaussSeidelForwardIteration(IVectorView rhsVector, IVector lhsVector);
	}
}
