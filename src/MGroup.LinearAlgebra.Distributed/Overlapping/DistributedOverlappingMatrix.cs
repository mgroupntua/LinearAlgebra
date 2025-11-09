using System;
using System.Collections.Concurrent;
using System.Collections.Generic;

using MGroup.Environments;
using MGroup.LinearAlgebra.Commons;
using MGroup.LinearAlgebra.Exceptions;
using MGroup.LinearAlgebra.Implementations;
using MGroup.LinearAlgebra.Matrices;
using MGroup.LinearAlgebra.Reduction;
using MGroup.LinearAlgebra.Vectors;

namespace MGroup.LinearAlgebra.Distributed.Overlapping
{
	public sealed class DistributedOverlappingMatrix<TMatrix> : DefaultMatrix
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

		public override int NumRows => Indexer.NumGlobalIndices;

		public override int NumColumns => Indexer.NumGlobalIndices;

		public override double this[int rowIdx, int colIdx]
		{
			get
			{
				Indexer.CheckGlobalIndex2D(rowIdx, colIdx);
				IReadOnlyDictionary<int, int> localRowIndices = Indexer.FindLocalIndicesOf(rowIdx);
				IReadOnlyDictionary<int, int> localColIndices = Indexer.FindLocalIndicesOf(colIdx);

				foreach ((int nodeID, int localColIdx) in localColIndices)
				{
					if (!LocalMatrices.ContainsKey(nodeID))
					{
						throw new Exception("This should not have happened. The distributed matrix is not created correctly.");
					}

					if (localRowIndices.TryGetValue(nodeID, out int localRowIdx))
					{
						return LocalMatrices[nodeID][localRowIdx, localColIdx]; // Return the first to be found.
					}
				}

				// If we reached this line, then the entry (rowIdx, colIdx) is not explicitly stored (structural zero)
				return 0.0;
			}
			set
			{
				Indexer.CheckGlobalIndex2D(rowIdx, colIdx);
				IReadOnlyDictionary<int, int> localRowIndices = Indexer.FindLocalIndicesOf(rowIdx);
				IReadOnlyDictionary<int, int> localColIndices = Indexer.FindLocalIndicesOf(colIdx);

				foreach ((int nodeID, int localColIdx) in localColIndices)
				{
					if (!LocalMatrices.ContainsKey(nodeID))
					{
						throw new Exception("This should not have happened. The distributed matrix is not created correctly.");
					}

					if (localRowIndices.TryGetValue(nodeID, out int localRowIdx))
					{
						LocalMatrices[nodeID].Set(localRowIdx, localColIdx, value); // Do this in all instances.
					}
				}

				// If we reached this line, then the entry (rowIdx, colIdx) is not explicitly stored (structural zero)
				if (value != 0.0)
				{
					throw new SparsityPatternModifiedException(
						$"The entry ({rowIdx}, {colIdx}) is a structural zero and cannot be changed");
				}
			}
		}

		public override void AxpyIntoThis(IMatrixView otherMatrix, double otherCoefficient)
		{
			if (otherMatrix is DistributedOverlappingMatrix<TMatrix> casted)
			{
				AxpyIntoThis(casted, otherCoefficient);
			}
			else
			{
				base.AxpyIntoThis(otherMatrix, otherCoefficient);
			}
		}

		public void AxpyIntoThis(DistributedOverlappingMatrix<TMatrix> otherMatrix, double otherCoefficient)
		{
			if (Indexer.IsCompatibleWith(otherMatrix.Indexer))
			{
				Environment.DoPerNode(
					node => this.LocalMatrices[node].AxpyIntoThis(otherMatrix.LocalMatrices[node], otherCoefficient));
			}
			else
			{
				base.AxpyIntoThis(otherMatrix, otherCoefficient);
			}
		}

		public override void Clear()
		{
			Environment.DoPerNode(nodeID => this.LocalMatrices[nodeID].Clear());
		}

		public override IMatrix Copy(bool copyIndexingData) => CopyAsDistributed(copyIndexingData);

		public DistributedOverlappingMatrix<TMatrix> CopyAsDistributed(bool copyIndexingData = false)
		{
			var indexerCloned = copyIndexingData ? Indexer.DeepCopy() : Indexer;
			var copy = new DistributedOverlappingMatrix<TMatrix>(indexerCloned);
			Environment.DoPerNode(nodeID => copy.LocalMatrices[nodeID] = (TMatrix)this.LocalMatrices[nodeID].Copy());
			return copy;
		}

		public override Matrix CopyToFullMatrix()
		{
			var result = Matrix.CreateZero(NumRows, NumColumns);
			Environment.DoPerNodeSerially(nodeID =>
			{
				TMatrix localMatrix = LocalMatrices[nodeID];
				for (int localRowIdx = 0; localRowIdx < localMatrix.NumRows; localRowIdx++)
				{
					int globalRowIdx = Indexer.FindGlobalIndexOf(nodeID, localRowIdx);
					for (int localColIdx = 0; localColIdx < localMatrix.NumColumns; localColIdx++)
					{
						int globalColIdx = Indexer.FindGlobalIndexOf(nodeID, localColIdx);
						result[globalRowIdx, globalColIdx] = localMatrix[localRowIdx, localColIdx];
					}
				}
			});
			return result;
		}

		public override IMatrix CreateZeroMatrixWithSameFormat() => CreateZeroMatrixSame();

		public DistributedOverlappingMatrix<TMatrix> CreateZeroMatrixSame()
			=> new DistributedOverlappingMatrix<TMatrix>(Indexer);

		public override void DoEntrywiseIntoThis(IMatrixView otherMatrix, Func<double, double, double> binaryOperation)
		{
			if (otherMatrix is DistributedOverlappingMatrix<TMatrix> casted)
			{
				DoEntrywiseIntoThis(casted, binaryOperation);
			}
			else
			{
				base.DoEntrywiseIntoThis(otherMatrix, binaryOperation);
			}
		}

		public void DoEntrywiseIntoThis(
			DistributedOverlappingMatrix<TMatrix> otherMatrix, Func<double, double, double> binaryOperation)
		{
			if (Indexer.IsCompatibleWith(otherMatrix.Indexer))
			{
				Environment.DoPerNode(
					node => this.LocalMatrices[node].DoEntrywiseIntoThis(otherMatrix.LocalMatrices[node], binaryOperation));
			}
			else
			{
				base.DoEntrywiseIntoThis(otherMatrix, binaryOperation);
			}
		}

		public override IMatrix DoToAllEntries(Func<double, double> unaryOperation)
		{
			DistributedOverlappingMatrix<TMatrix> result = CopyAsDistributed();
			result.DoToAllEntriesIntoThis(unaryOperation);
			return result;
		}

		public override void DoToAllEntriesIntoThis(Func<double, double> unaryOperation)
		{
			Environment.DoPerNode(node => this.LocalMatrices[node].DoToAllEntriesIntoThis(unaryOperation));
		}

		public override bool Equals(IIndexable2D other, double tolerance = 1E-7)
		{
			if (other is DistributedOverlappingMatrix<TMatrix> casted)
			{
				return this.Equals(casted, tolerance);
			}

			return false;
		}

		public bool Equals(DistributedOverlappingMatrix<TMatrix> other, double tolerance = 1E-7)
		{
			if (!this.Indexer.IsCompatibleWith(other.Indexer))
			{
				return false;
			}

			Dictionary<int, bool> flags = Environment.CalcNodeData(
					node => this.LocalMatrices[node].Equals(other.LocalMatrices[node], tolerance));
			return Environment.AllReduceAnd(flags);
		}

		public override bool HasSameFormat(IMatrixView otherMatrix)
		{
			if (otherMatrix is DistributedOverlappingMatrix<TMatrix> casted)
			{
				return this.Indexer.IsCompatibleWith(casted.Indexer);
			}

			return false;
		}

		public bool HasSameFormat(DistributedOverlappingMatrix<TMatrix> otherMatrix)
			=> this.Indexer.IsCompatibleWith(otherMatrix.Indexer);

		public override void LinearCombinationIntoThis(double thisCoefficient, IMatrixView otherMatrix, double otherCoefficient)
		{
			if (otherMatrix is DistributedOverlappingMatrix<TMatrix> casted)
			{
				LinearCombinationIntoThis(thisCoefficient, casted, otherCoefficient);
			}
			else
			{
				base.LinearCombinationIntoThis(thisCoefficient, otherMatrix, otherCoefficient);
			}
		}

		public void LinearCombinationIntoThis(
			double thisCoefficient, DistributedOverlappingMatrix<TMatrix> otherMatrix, double otherCoefficient)
		{
			if (Indexer.IsCompatibleWith(otherMatrix.Indexer))
			{
				Environment.DoPerNode(
					node => this.LocalMatrices[node].LinearCombinationIntoThis(
						thisCoefficient, otherMatrix.LocalMatrices[node], otherCoefficient));
			}
			else
			{
				base.LinearCombinationIntoThis(thisCoefficient, otherMatrix, otherCoefficient);
			}
		}

		public override IVector Multiply(IVectorView vector, bool transposeThis = false)
		{
			if (vector is DistributedOverlappingVector lhsCasted)
			{
				DistributedOverlappingVector result = lhsCasted.CreateZeroVectorSame();
				MultiplyIntoResult(lhsCasted, result, transposeThis);
				return result;
			}
			else
			{
				return base.Multiply(vector, transposeThis);
			}
		}

		public override void MultiplyIntoResult(IVectorView lhsVector, IVector rhsVector, bool transposeThis = false)
		{
			if ((lhsVector is DistributedOverlappingVector lhsCasted) && (rhsVector is DistributedOverlappingVector rhsCasted))
			{
				MultiplyIntoResult(lhsCasted, rhsCasted, transposeThis);
			}
			else
			{
				base.MultiplyIntoResult(lhsVector, rhsVector, transposeThis);
			}
		}

		public void MultiplyIntoResult(DistributedOverlappingVector lhsVector, DistributedOverlappingVector rhsVector,
			bool transposeThis = false)
		{
			if (this.Indexer.IsCompatibleWith(lhsVector.Indexer) && this.Indexer.IsCompatibleWith(rhsVector.Indexer))
			{
				Action<int> multiplyLocal = nodeID =>
				{
					TMatrix localA = this.LocalMatrices[nodeID];
					Vector localX = lhsVector.LocalVectors[nodeID];
					Vector localY = rhsVector.LocalVectors[nodeID];
					localA.MultiplyIntoResult(localX, localY, transposeThis);
				};
				Environment.DoPerNode(multiplyLocal);

				rhsVector.SumOverlappingEntries();
			}
			else
			{
				base.MultiplyIntoResult(lhsVector, rhsVector, transposeThis);
			}
		}

		public override IMatrix Scale(double scalar)
		{
			DistributedOverlappingMatrix<TMatrix> result = CopyAsDistributed();
			result.ScaleIntoThis(scalar);
			return result;
		}

		public override void ScaleIntoThis(double coefficient)
		{
			Environment.DoPerNode(nodeID => this.LocalMatrices[nodeID].ScaleIntoThis(coefficient));
		}

		public override IMatrix Transpose() => TransposeDistributed();

		public DistributedOverlappingMatrix<TMatrix> TransposeDistributed()
		{
			var transpose = new DistributedOverlappingMatrix<TMatrix>(Indexer);
			Environment.DoPerNode(nodeID => transpose.LocalMatrices[nodeID] = (TMatrix)this.LocalMatrices[nodeID].Transpose());
			return transpose;
		}
	}
}
