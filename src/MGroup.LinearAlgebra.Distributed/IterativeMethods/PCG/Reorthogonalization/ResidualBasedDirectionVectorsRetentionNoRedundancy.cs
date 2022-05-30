using System;
using System.Collections.Generic;
using System.Text;

namespace MGroup.LinearAlgebra.Distributed.IterativeMethods.PCG.Reorthogonalization
{
	public class ResidualBasedDirectionVectorsRetentionNoRedundancy : IDirectionVectorsRetention
	{
		private readonly double minResidualNormRatioToKeep;

		private ReorthogonalizedPcg pcg;

		public ResidualBasedDirectionVectorsRetentionNoRedundancy(double minResidualNormRatioToKeep)
		{
			this.minResidualNormRatioToKeep = minResidualNormRatioToKeep;
		}

		public void DiscardDirectionVectors() { }

		public void Intialize(ReorthogonalizedPcg pcg)
		{
			this.pcg = pcg;
		}

		public bool KeepUsingReorthogonalization()
		{
			if (double.IsNaN(pcg.ResidualNormRatio))
			{
				return true; // First iteration - not initialized yet.
			}
			else
			{
				return pcg.ResidualNormRatio > minResidualNormRatioToKeep;
			}
		}
	}
}
