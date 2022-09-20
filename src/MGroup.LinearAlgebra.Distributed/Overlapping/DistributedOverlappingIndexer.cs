using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;
using MGroup.Environments;
using MGroup.LinearAlgebra.Distributed.LinearAlgebraExtensions;
using MGroup.LinearAlgebra.Matrices;
using MGroup.Solvers.DDM.LinearSystem;
using MGroup.MSolve.Solution.LinearSystem;

namespace MGroup.LinearAlgebra.Distributed.Overlapping
{
	/// <summary>
	/// Manages the indices for a <see cref="DistributedOverlappingVector"/>, <see cref="DistributedOverlappingVector"/>, etc.
	/// Supports multiple local vectors, each of which may have none, some or all its entries in common with other local vectors.
	/// Specifies the relationships between these common entries. When dealing with multiple distributed vectors that have the 
	/// same indexing pattern, reuse the same instance of <see cref="DistributedOverlappingIndexer"/>.
	/// </summary>
	/// <remarks>
	/// In interface problems of PSM and related DDMs, local vectors have all their entries in common with other local vectors, 
	/// since only boundary dofs take part in the interface problem. In GSI or a GSI-like treatment of other DDMs' coarse 
	/// problems, local vectors have some of their entries (boundary dofs) in common with other local vectors, while the rest
	/// entries (internal dofs) are unique for each local vector. 
	/// </remarks>
	public class DistributedOverlappingIndexer : IDistributedIndexer
	{
		private readonly Dictionary<int, Local> localIndexers;
		private int numUniqueEntries = int.MinValue;

		public DistributedOverlappingIndexer(IComputeEnvironment environment)
		{
			localIndexers = environment.CalcNodeData(
				n => new Local(environment.GetComputeNode(n)));
			Environment = environment;
		}

		public IComputeEnvironment Environment { get; }

		public int CountUniqueEntries()
		{
			if (numUniqueEntries == int.MinValue)
			{
				Dictionary<int, double> countPerNode = Environment.CalcNodeData(node =>
				{
					double[] inverseMultiplicities = localIndexers[node].InverseMultiplicities;

					double localCount = 0.0;
					for (int i = 0; i < inverseMultiplicities.Length; ++i)
					{
						localCount += inverseMultiplicities[i];
					}

					return localCount;
				});
				double globalCount = Environment.AllReduceSum(countPerNode);
				numUniqueEntries = (int)globalCount;
			}
			return numUniqueEntries;
		}

		public DistributedOverlappingIndexer.Local GetLocalComponent(int nodeID) => localIndexers[nodeID];

		public bool IsCompatibleWith(IDistributedIndexer other) => this == other;

		public DistributedOverlappingMatrix<TMatrix> CheckCompatibleMatrix<TMatrix>(IGlobalMatrix matrix)
			where TMatrix : class, IMatrix
		{
			if (matrix is DistributedOverlappingMatrix<TMatrix> distributedMatrix)
			{
				if (matrix.CheckForCompatibility == false || distributedMatrix.Indexer == this)
				{
					return distributedMatrix;
				}
			}

			throw new NonMatchingFormatException("The provided matrix has a different format than this indexer. " +
				"Their entries correspond to different dofs or they are distributed differently across compute nodes");
		}

		public DistributedOverlappingVector CheckCompatibleVector(IGlobalVector vector)
		{
			if (vector is DistributedOverlappingVector distributedVector)
			{
				if (vector.CheckForCompatibility == false || distributedVector.Indexer == this)
				{
					return distributedVector;
				}
			}

			//TODO: Perhaps the error msg or even the whole check can be injected into the constructor.
			throw new NonMatchingFormatException("The provided vector has a different format than this indexer. " +
				"Their entries correspond to different dofs or they are distributed differently across compute nodes");
		}

		public bool IsCompatibleVector(DistributedOverlappingVector vector) => vector.Indexer == this;

		/// <summary>
		/// All indexing data and functionality of <see cref="DistributedOverlappingIndexer"/>, but only for the local vector, 
		/// matrix, etc. that corresponds to a specific <see cref="ComputeNode"/>.
		/// </summary>
		public class Local
		{
			private Dictionary<int, int[]> commonEntriesWithNeighbors;

			public Local(ComputeNode node)
			{
				this.Node = node;
			}

			/// <summary>
			/// Neighboring <see cref="ComputeNode"/>s of this <see cref="Node"/> with local vectors that have at least 1 common 
			/// entry with the local vector of this <see cref="Node"/>.
			/// </summary>
			public SortedSet<int> ActiveNeighborsOfNode { get; private set; }

			public double[] InverseMultiplicities { get; private set; } 

			public ComputeNode Node { get; }

			public int NumEntries { get; private set; }

			//TODO: cache a buffer for sending and a buffer for receiving inside Indexer (lazily or not) and just return them. 
			//      Also provide an option to request newly initialized buffers. It may be better to have dedicated Buffer classes to
			//      handle all that logic (e.g. keeping allocated buffers in a LinkedList, giving them out & locking them, 
			//      freeing them in clients, etc.
			public ConcurrentDictionary<int, double[]> CreateBuffersForAllToAllWithNeighbors()
			{
				//TODOMPI: dictionaries that contain per node values should be requested from the environment, which knows their
				//      type (Dictionary/ConcurrentDictionary), capacity and concurrency level.
				var buffers = new ConcurrentDictionary<int, double[]>(); 
				foreach (int neighborID in ActiveNeighborsOfNode)
				{
					buffers[neighborID] = new double[commonEntriesWithNeighbors[neighborID].Length];
				}
				return buffers;
			}

			public void Initialize(int numTotalEntries, Dictionary<int, int[]> commonEntriesWithNeighbors)
			{
				this.NumEntries = numTotalEntries;
				ActiveNeighborsOfNode = new SortedSet<int>(commonEntriesWithNeighbors.Keys);
				Debug.Assert(Node.Neighbors.IsSupersetOf(ActiveNeighborsOfNode));
				this.commonEntriesWithNeighbors = commonEntriesWithNeighbors;
				FindMultiplicities();
			}

			/// <summary>
			/// Copy data shallowly from <paramref name="other"/>.
			/// </summary>
			/// <param name="other"></param>
			public void InitializeFrom(Local other)
			{
				this.NumEntries = other.NumEntries;
				this.commonEntriesWithNeighbors = other.commonEntriesWithNeighbors;
				this.ActiveNeighborsOfNode = other.ActiveNeighborsOfNode;
				this.InverseMultiplicities = other.InverseMultiplicities;
			}

			public int[] GetCommonEntriesWithNeighbor(int neighbor) => commonEntriesWithNeighbors[neighbor];

			public void FindMultiplicities()
			{
				var multiplicities = new int[NumEntries];
				for (int i = 0; i < NumEntries; ++i) multiplicities[i] = 1;
				foreach (int[] commonEntries in commonEntriesWithNeighbors.Values)
				{
					foreach (int i in commonEntries) multiplicities[i] += 1;
				}

				InverseMultiplicities = new double[NumEntries];
				for (int i = 0; i < NumEntries; ++i) InverseMultiplicities[i] = 1.0 / multiplicities[i];
			}
		}
	}
}
