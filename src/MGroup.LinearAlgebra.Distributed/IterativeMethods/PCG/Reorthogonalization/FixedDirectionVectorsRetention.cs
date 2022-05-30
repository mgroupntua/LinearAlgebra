using System;
using System.Collections.Generic;
using System.Text;

namespace MGroup.LinearAlgebra.Distributed.IterativeMethods.PCG.Reorthogonalization
{
	public class FixedDirectionVectorsRetention : IDirectionVectorsRetention
	{
		private readonly int numVectorsToKeep;
		private readonly bool keepFirstVectors;

		private ReorthogonalizedPcg pcg;

		public FixedDirectionVectorsRetention(int numVectorsToKeep, bool keepOldestVectors = true)
		{
			this.numVectorsToKeep = numVectorsToKeep;
			this.keepFirstVectors = keepOldestVectors;
		}

		public void DiscardDirectionVectors()
		{
			int numVectorsToDiscard = pcg.ReorthoCache.Directions.Count - numVectorsToKeep;
			if (keepFirstVectors)
			{
				pcg.ReorthoCache.RemoveNewDirectionVectorData(numVectorsToDiscard);
			}
			else
			{
				pcg.ReorthoCache.RemoveOldDirectionVectorData(numVectorsToDiscard);
			}
		}

		public void Intialize(ReorthogonalizedPcg pcg)
		{
			this.pcg = pcg;
		}

		public bool KeepUsingReorthogonalization() => true;
	}
}
