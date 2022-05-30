using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;
using MGroup.LinearAlgebra.Vectors;
using MGroup.Environments;
using MGroup.LinearAlgebra.Distributed.Overlapping;
using MGroup.MSolve.Solution.LinearSystem;

//TODOMPI: Needs testing
namespace MGroup.LinearAlgebra.Distributed.IterativeMethods.Preconditioning
{
	public class DistributedOverlappingJacobiPreconditioner : IPreconditioner
	{
		private readonly IComputeEnvironment environment;

		/// <summary>
		/// Creates a Jacobi preconditioner in distributed environment.
		/// </summary>
		/// <param name="environment">
		/// The computing environment that will be used for the operations during this constructor and during 
		/// <see cref="Apply(IGlobalVector, IGlobalVector)"/>.
		/// </param>
		/// <param name="diagonal">
		/// A distributed vector that contains the diagonal entries of each local matrix that corresponds to a 
		/// <see cref="ComputeNode"/> of the <paramref name="environment"/>. If an entry is overlapping, namely if it exists
		/// in many neighboring local diagonal vectors, then its value must be the same in all these local vectors.
		/// </param>
		public DistributedOverlappingJacobiPreconditioner(IComputeEnvironment environment,
			DistributedOverlappingVector diagonal)
		{
			this.environment = environment;
			this.Indexer = diagonal.Indexer;

			Func<int, Vector> invertDiagonal = nodeID => diagonal.LocalVectors[nodeID].DoToAllEntries(x => 1 / x);
			this.LocalInverseDiagonals = environment.CalcNodeData(invertDiagonal);
		}

		public IDistributedIndexer Indexer { get; }

		public Dictionary<int, Vector> LocalInverseDiagonals { get; }

		public void Apply(IGlobalVector input, IGlobalVector output)
		{
			if ((input is DistributedOverlappingVector lhsCasted) && (output is DistributedOverlappingVector rhsCasted))
			{
				Multiply(lhsCasted, rhsCasted);
			}
			else
			{
				throw new ArgumentException(
					"This operation is legal only if the left-hand-side and righ-hand-side vectors are distributed" +
					" with overlapping entries.");
			}
		}

		public void Multiply(DistributedOverlappingVector input, DistributedOverlappingVector output)
		{
			//TODOMPI: also check that environment is the same between M,x and M,y
			Debug.Assert(this.Indexer.IsCompatibleWith(input.Indexer) /*&& (this.environment == lhs.environment)*/);
			Debug.Assert(this.Indexer.IsCompatibleWith(output.Indexer) /*&& (this.environment == rhs.environment)*/);

			Action<int> multiplyLocal = nodeID =>
			{
				Vector localX = input.LocalVectors[nodeID];
				Vector localY = output.LocalVectors[nodeID];
				Vector localDiagonal = this.LocalInverseDiagonals[nodeID];
				localY.CopyFrom(localX);
				localY.MultiplyEntrywiseIntoThis(localDiagonal);
			};
			environment.DoPerNode(multiplyLocal);

			//TODOMPI: do we need to call output.SumOverlappingEntries() here? Is this need covered by the fact that 
			//      LocalInverseDiagonals already have the total stiffnesses?
		}
	}
}
