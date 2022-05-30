using System;
using System.Collections.Generic;
using System.Text;

namespace MGroup.LinearAlgebra.Distributed.IterativeMethods.PCG.Reorthogonalization
{
	public interface IDirectionVectorsRetention
	{
		void DiscardDirectionVectors();

		void Intialize(ReorthogonalizedPcg pcg);

		bool KeepUsingReorthogonalization();
	}
}
