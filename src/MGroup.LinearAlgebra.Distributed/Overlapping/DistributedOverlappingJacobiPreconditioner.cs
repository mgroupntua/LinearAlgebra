//TODOMPI: Needs testing
namespace MGroup.LinearAlgebra.Distributed.Overlapping
{
	using System;
	using System.Collections.Generic;
	using System.Diagnostics;
	using System.Text;
	using MGroup.LinearAlgebra.Vectors;
	using MGroup.Environments;
	using MGroup.LinearAlgebra.Matrices;
	using MGroup.LinearAlgebra.Iterative.Preconditioning;

	public class DistributedOverlappingJacobiPreconditioner : IPreconditioner
	{
		private readonly IComputeEnvironment environment;
		private readonly DistributedOverlappingVector diagonal;

		/// <summary>
		/// Creates a Jacobi preconditioner in distributed environment.
		/// </summary>
		/// <param name="environment">
		/// The computing environment that will be used for the operations during this constructor and during 
		/// <see cref="Apply(IVector, IVector)"/>.
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
			this.diagonal = diagonal;
			Indexer = diagonal.Indexer;

			Func<int, Vector> invertDiagonal = nodeID => diagonal.LocalVectors[nodeID].DoToAllEntries(x => 1 / x);
			LocalInverseDiagonals = environment.CalcNodeData(invertDiagonal);
		}

		public IDistributedIndexer Indexer { get; }

		public Dictionary<int, Vector> LocalInverseDiagonals { get; }

		public void Apply(DistributedOverlappingVector input, DistributedOverlappingVector output)
		{
			//TODOMPI: also check that environment is the same between M,x and M,y
			Debug.Assert(Indexer.IsCompatibleWith(input.Indexer) /*&& (this.environment == lhs.environment)*/);
			Debug.Assert(Indexer.IsCompatibleWith(output.Indexer) /*&& (this.environment == rhs.environment)*/);

			Action<int> multiplyLocal = nodeID =>
			{
				var localX = input.LocalVectors[nodeID];
				var localY = output.LocalVectors[nodeID];
				var localDiagonal = LocalInverseDiagonals[nodeID];
				localY.CopyFrom(localX);
				localY.MultiplyEntrywiseIntoThis(localDiagonal);
			};
			environment.DoPerNode(multiplyLocal);

			//TODOMPI: do we need to call output.SumOverlappingEntries() here? Is this need covered by the fact that 
			//      LocalInverseDiagonals already have the total stiffnesses?
		}

		public IPreconditioner CopyWithInitialSettings() => new DistributedOverlappingJacobiPreconditioner(environment, diagonal);


		public void SolveLinearSystem(IReadOnlyVector rhsVector, IVector lhsVector)
		{
			if (rhsVector is DistributedOverlappingVector rhsCasted && lhsVector is DistributedOverlappingVector lhsCasted)
			{
				Apply(rhsCasted, lhsCasted);
			}
			else
			{
				throw new ArgumentException(
					"This operation is legal only if the left-hand-side and righ-hand-side vectors are distributed" +
					" with overlapping entries.");
			}
		}

		public void UpdateMatrix(IReadOnlyMatrix matrix, bool isPatternModified) { } 
	}
}
