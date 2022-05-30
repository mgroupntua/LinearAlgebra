using System;
using System.Collections.Generic;
using System.Text;

namespace MGroup.LinearAlgebra.Distributed.IterativeMethods.PCG.Reorthogonalization
{
	public class PercentageDirectionVectorsRetention : IDirectionVectorsRetention
	{
		private readonly double percentageOfVectorsToKeep;
		private readonly bool keepFirstVectors;

		private ReorthogonalizedPcg pcg;

		public PercentageDirectionVectorsRetention(double percentageOfVectorsToKeep, bool keepOldestVectors = true)
		{
			this.percentageOfVectorsToKeep = percentageOfVectorsToKeep;
			this.keepFirstVectors = keepOldestVectors;
		}
		public void DiscardDirectionVectors()
		{
			int numVectorsTotal = pcg.ReorthoCache.Directions.Count;
			int numVectorsToKeep = (int)Math.Round(percentageOfVectorsToKeep * numVectorsTotal);
			if (numVectorsToKeep <= 0)
			{
				pcg.ReorthoCache.Clear();
			}
			else if (numVectorsToKeep >= numVectorsTotal)
			{
				return;
			}
			else
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
		}

		public void Intialize(ReorthogonalizedPcg pcg)
		{
			this.pcg = pcg;
		}

		public bool KeepUsingReorthogonalization() => true;
	}
}
