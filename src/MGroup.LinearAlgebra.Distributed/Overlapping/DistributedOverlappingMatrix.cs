using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;
using MGroup.Environments;
using MGroup.MSolve.Solution.LinearSystem;
using MGroup.LinearAlgebra.Distributed.LinearAlgebraExtensions;
using MGroup.LinearAlgebra.Distributed.Overlapping;
using MGroup.LinearAlgebra.Matrices;
using MGroup.LinearAlgebra.Vectors;
using MGroup.MSolve.Solution.LinearSystem;

namespace MGroup.Solvers.DDM.LinearSystem
{
	public class DistributedOverlappingMatrix<TMatrix> : IGlobalMatrix
		where TMatrix : class, IMatrix
	{
		public DistributedOverlappingMatrix(DistributedOverlappingIndexer indexer)
		{
			this.Indexer = indexer;
			this.Environment = indexer.Environment;
		}

		public IComputeEnvironment Environment { get; }

		public DistributedOverlappingIndexer Indexer { get; }

		public ConcurrentDictionary<int, TMatrix> LocalMatrices { get; } = new ConcurrentDictionary<int, TMatrix>();
		public bool CheckForCompatibility { get; set; } = true;

		public void AxpyIntoThis(IGlobalMatrix otherMatrix, double otherCoefficient)
		{
			DistributedOverlappingMatrix<TMatrix> distributedOther = Indexer.CheckCompatibleMatrix<TMatrix>(otherMatrix);
			Action<int> localOperation = nodeID =>
			{
				TMatrix thisSubdomainMatrix = this.LocalMatrices[nodeID];
				TMatrix otherSubdomainMatrix = distributedOther.LocalMatrices[nodeID];
				thisSubdomainMatrix.AxpyIntoThis(otherSubdomainMatrix, otherCoefficient);
			};
			Environment.DoPerNode(localOperation);
		}

		public void Clear()
		{
			Environment.DoPerNode(nodeID => this.LocalMatrices[nodeID].Clear());
		}

		public IGlobalMatrix Copy()
		{
			var copy = new DistributedOverlappingMatrix<TMatrix>(Indexer);
			Environment.DoPerNode(nodeID => copy.LocalMatrices[nodeID] = (TMatrix)this.LocalMatrices[nodeID].Copy());
			return copy;
		}

		public void LinearCombinationIntoThis(double thisCoefficient, IGlobalMatrix otherMatrix, double otherCoefficient)
		{
			DistributedOverlappingMatrix<TMatrix> distributedOther = Indexer.CheckCompatibleMatrix<TMatrix>(otherMatrix);
			Action<int> localOperation = nodeID =>
			{
				TMatrix thisSubdomainMatrix = this.LocalMatrices[nodeID];
				TMatrix otherSubdomainMatrix = distributedOther.LocalMatrices[nodeID];
				thisSubdomainMatrix.LinearCombinationIntoThis(thisCoefficient, otherSubdomainMatrix, otherCoefficient);
			};
			Environment.DoPerNode(localOperation);
		}

		public void MultiplyVector(IGlobalVector input, IGlobalVector output)
		{
			DistributedOverlappingVector distributedInput = Indexer.CheckCompatibleVector(input);
			DistributedOverlappingVector distributedOutput = Indexer.CheckCompatibleVector(output);
			MultiplyVector(distributedInput, distributedOutput);
		}

		public void MultiplyVector(DistributedOverlappingVector input, DistributedOverlappingVector output)
		{
			Indexer.CheckCompatibleVector(input);
			Indexer.CheckCompatibleVector(output);

			Action<int> multiplyLocal = nodeID =>
			{
				TMatrix localA = this.LocalMatrices[nodeID];
				Vector localX = input.LocalVectors[nodeID];
				Vector localY = output.LocalVectors[nodeID];
				localA.MultiplyIntoResult(localX, localY);
			};
			Environment.DoPerNode(multiplyLocal);

			output.SumOverlappingEntries();
		}

		public void ScaleIntoThis(double coefficient)
		{
			Environment.DoPerNode(nodeID => this.LocalMatrices[nodeID].ScaleIntoThis(coefficient));
		}
	}
}
