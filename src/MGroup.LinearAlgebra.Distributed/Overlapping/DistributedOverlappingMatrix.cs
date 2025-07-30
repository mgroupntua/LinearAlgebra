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

using static MGroup.LinearAlgebra.Distributed.Overlapping.CompatibilityUtilities;

namespace MGroup.LinearAlgebra.Distributed.Overlapping
{
	public class DistributedOverlappingMatrix<TMatrix> : IMatrix //IGlobalMatrix
		where TMatrix : class, IMatrix
	{
		private GlobalIndexer globalIndexer;

		public DistributedOverlappingMatrix(DistributedOverlappingIndexer indexer)
		{
			this.Indexer = indexer;
			this.Environment = indexer.Environment;
		}

		public IComputeEnvironment Environment { get; }

		public DistributedOverlappingIndexer Indexer { get; }

		public ConcurrentDictionary<int, TMatrix> LocalMatrices { get; } = new ConcurrentDictionary<int, TMatrix>();

		public int NumRows => Indexer.NumUniqueEntries;

		public int NumColumns => Indexer.NumUniqueEntries;

		public MatrixSymmetry MatrixSymmetry { get; set; } = MatrixSymmetry.Unknown;

		public double this[int rowIdx, int colIdx]
		{
			get
			{
				CreateGlobalIndexerIfMissing();
				globalIndexer.CheckGlobalIndex2D(rowIdx, colIdx);
				IReadOnlyDictionary<int, int> localRowIndices = globalIndexer.FindLocalIndicesOf(rowIdx);
				IReadOnlyDictionary<int, int> localColIndices =  globalIndexer.FindLocalIndicesOf(colIdx);

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
		}

		public void AxpyIntoThis(IMatrixView otherMatrix, double otherCoefficient)
		{
			DistributedOverlappingMatrix<TMatrix> distributedOther = CastToDistributed<TMatrix>(otherMatrix);
			CheckSameFormat(this,distributedOther);
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


		IMatrix IMatrixView.Copy(bool copyIndexingData) => Copy(copyIndexingData);

		public DistributedOverlappingMatrix<TMatrix> Copy(bool copyIndexingData = false)
		{
			var indexerCloned = copyIndexingData ? Indexer.DeepCopy() : Indexer;
			var copy = new DistributedOverlappingMatrix<TMatrix>(indexerCloned);
			Environment.DoPerNode(nodeID => copy.LocalMatrices[nodeID] = (TMatrix)this.LocalMatrices[nodeID].Copy());
			return copy;
		}

		public Matrix CopyToFullMatrix()
		{
			CreateGlobalIndexerIfMissing();
			var result = Matrix.CreateZero(NumRows, NumColumns);
			Environment.DoPerNodeSerially(nodeID =>
			{
				TMatrix localMatrix = LocalMatrices[nodeID];
				for (int localRowIdx = 0; localRowIdx < localMatrix.NumRows; localRowIdx++)
				{
					int globalRowIdx = globalIndexer.FindGlobalIndexOf(nodeID, localRowIdx);
					for (int localColIdx = 0; localColIdx < localMatrix.NumColumns; localColIdx++)
					{
						int globalColIdx = globalIndexer.FindGlobalIndexOf(nodeID, localColIdx);
						result[globalRowIdx, globalColIdx] = localMatrix[localRowIdx, localColIdx];
					}
				}
			});
			return result;
		}

		public DistributedOverlappingMatrix<TMatrix> CreateZeroMatrixWithSameFormat() //TODO: Move this to IMatrixView
			=> new DistributedOverlappingMatrix<TMatrix>(Indexer);

		public IMatrix DoEntrywise(IMatrixView other, Func<double, double, double> binaryOperation)
		{
			DistributedOverlappingMatrix<TMatrix> result = Copy();
			result.DoEntrywiseIntoThis(other, binaryOperation);
			return result;
		}

		public void DoEntrywiseIntoThis(IMatrixView other, Func<double, double, double> binaryOperation)
			=> DoEntrywiseIntoThis(CastToDistributed<TMatrix>(other), binaryOperation);

		public void DoEntrywiseIntoThis(DistributedOverlappingMatrix<TMatrix> other, Func<double, double, double> binaryOperation)
		{
			CheckSameFormat(this, other);
			Environment.DoPerNode(node => this.LocalMatrices[node].DoEntrywiseIntoThis(other, binaryOperation));
		}

		public IMatrix DoToAllEntries(Func<double, double> unaryOperation)
		{
			DistributedOverlappingMatrix<TMatrix> result = Copy();
			result.DoToAllEntriesIntoThis(unaryOperation);
			return result;
		}

		public void DoToAllEntriesIntoThis(Func<double, double> unaryOperation)
		{
			Environment.DoPerNode(node => this.LocalMatrices[node].DoToAllEntriesIntoThis(unaryOperation));
		}

		public bool Equals(IIndexable2D other, double tolerance = 1E-7)
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

		public Vector GetColumn(int colIndex) //TODO: IVectorView should implement ISliceable1D. Use default interface implementations
		{
			Preconditions.CheckIndexCol(this, colIndex);
			double[] columnVector = new double[NumRows];
			for (int i = 0; i < NumRows; i++)
			{
				columnVector[i] = this[i, colIndex];
			}

			return Vector.CreateFromArray(columnVector, false);
		}

		public Vector GetRow(int rowIndex)
		{
			Preconditions.CheckIndexRow(this, rowIndex);
			double[] rowVector = new double[NumColumns];
			for (int j = 0; j < NumColumns; j++)
			{
				rowVector[j] = this[rowIndex, j];
			}

			return Vector.CreateFromArray(rowVector, false);
		}

		public IMatrix GetSubmatrix(int[] rowIndices, int[] colIndices)
			=> DenseStrategies.GetSubmatrix(this, rowIndices, colIndices);
		
		public IMatrix GetSubmatrix(int rowStartInclusive, int rowEndExclusive, int colStartInclusive, int colEndExclusive)
			=> DenseStrategies.GetSubmatrix(this, rowStartInclusive, rowEndExclusive, colStartInclusive, colEndExclusive);

		public bool HasSameFormat(IIndexable2D other) //TODO: Move this to IMatrixView. Same for IVectorView
		{
			if (other is DistributedOverlappingVector casted)
			{
				return this.Indexer.IsCompatibleWith(casted.Indexer);
			}

			return false;
		}

		public void LinearCombinationIntoThis(double thisCoefficient, IMatrixView otherMatrix, double otherCoefficient)
		{
			DistributedOverlappingMatrix<TMatrix> distributedOther = CastToDistributed<TMatrix>(otherMatrix);
			CheckSameFormat(this, distributedOther);
			Action<int> localOperation = nodeID =>
			{
				TMatrix thisSubdomainMatrix = this.LocalMatrices[nodeID];
				TMatrix otherSubdomainMatrix = distributedOther.LocalMatrices[nodeID];
				thisSubdomainMatrix.LinearCombinationIntoThis(thisCoefficient, otherSubdomainMatrix, otherCoefficient);
			};
			Environment.DoPerNode(localOperation);
		}

		public Matrix MultiplyLeft(IMatrixView other, bool transposeThis = false, bool transposeOther = false) 
			=> other.MultiplyRight(this, transposeOther, transposeOther); //TODO: default interface implementation

		public Matrix MultiplyRight(IMatrixView other, bool transposeThis = false, bool transposeOther = false) //TODO: default interface implementation
		{
			if (transposeThis)
			{
				if (transposeOther)
				{
					Preconditions.CheckMultiplicationDimensions(this.NumRows, other.NumColumns);
					var result = Matrix.CreateZero(this.NumColumns, other.NumRows);
					for (int i = 0; i < this.NumColumns; i++)
					{
						for (int j = 0; j < other.NumRows; j++)
						{
							for (int k = 0; k < this.NumRows; k++)
							{
								result[i, j] += this[k, i] * other[j, k];
							}
						}
					}
					return result;
				}
				else
				{
					Preconditions.CheckMultiplicationDimensions(this.NumRows, other.NumRows);
					var result = Matrix.CreateZero(this.NumColumns, other.NumColumns);
					for (int i = 0; i < this.NumColumns; i++)
					{
						for (int j = 0; j < other.NumColumns; j++)
						{
							for (int k = 0; k < this.NumRows; k++)
							{
								result[i, j] += this[k, i] * other[k, j];
							}
						}
					}
					return result;
				}
			}
			else
			{
				if (transposeOther)
				{
					Preconditions.CheckMultiplicationDimensions(this.NumColumns, other.NumColumns);
					var result = Matrix.CreateZero(this.NumRows, other.NumRows);
					for (int i = 0; i < this.NumRows; i++)
					{
						for (int j = 0; j < other.NumRows; j++)
						{
							for (int k = 0; k < this.NumColumns; k++)
							{
								result[i, j] += this[i, k] * other[j, k];
							}
						}
					}
					return result;
				}
				else
				{
					Preconditions.CheckMultiplicationDimensions(this.NumColumns, other.NumRows);
					var result = Matrix.CreateZero(this.NumRows, other.NumColumns);
					for (int i = 0; i < this.NumRows; i++)
					{
						for (int j = 0; j < other.NumColumns; j++)
						{
							for (int k = 0; k < this.NumColumns; k++)
							{
								result[i, j] += this[i, k] * other[k, j];
							}
						}
					}
					return result;
				}
			}
		}

		public IVector Multiply(IVectorView vector, bool transposeThis = false) //TODO: Rename to MultiplyVector. Also IMatrixView must implement ILinearTransformation
		{
			DistributedOverlappingVector distributedLhs = CastToDistributed(vector);
			DistributedOverlappingVector result = distributedLhs.CreateZeroVectorWithSameFormat();
			MultiplyIntoResult(distributedLhs, result, transposeThis);
			return result;
		}

		public void MultiplyIntoResult(IVectorView lhsVector, IVector rhsVector, bool transposeThis = false)
		{
			DistributedOverlappingVector distributedLhs = CastToDistributed(lhsVector);
			DistributedOverlappingVector distributedRhs = CastToDistributed(rhsVector);
			MultiplyIntoResult(distributedLhs, distributedRhs, transposeThis);
		}

		public void MultiplyIntoResult(DistributedOverlappingVector lhsVector, DistributedOverlappingVector rhsVector,
			bool transposeThis = false)
		{
			if (transposeThis)
			{
				throw new NotImplementedException();
			}

			CheckSameFormat(this, lhsVector);
			CheckSameFormat(this, rhsVector);

			Action<int> multiplyLocal = nodeID =>
			{
				TMatrix localA = this.LocalMatrices[nodeID];
				Vector localX = lhsVector.LocalVectors[nodeID];
				Vector localY = rhsVector.LocalVectors[nodeID];
				localA.MultiplyIntoResult(localX, localY);
			};
			Environment.DoPerNode(multiplyLocal);

			rhsVector.SumOverlappingEntries();
		}

		public double Reduce(double identityValue, ProcessEntry processEntry, ProcessZeros processZeros, Finalize finalize)
			=> throw new NotImplementedException("Environment must define these reductions"); //TODO: use default implementation

		public void ScaleIntoThis(double coefficient)
		{
			Environment.DoPerNode(nodeID => this.LocalMatrices[nodeID].ScaleIntoThis(coefficient));
		}

		public void SetEntryRespectingPattern(int rowIdx, int colIdx, double value)
		{
			CreateGlobalIndexerIfMissing();
			globalIndexer.CheckGlobalIndex2D(rowIdx, colIdx);
			IReadOnlyDictionary<int, int> localRowIndices = globalIndexer.FindLocalIndicesOf(rowIdx);
			IReadOnlyDictionary<int, int> localColIndices = globalIndexer.FindLocalIndicesOf(colIdx);

			foreach ((int nodeID, int localColIdx) in localColIndices)
			{
				if (!LocalMatrices.ContainsKey(nodeID))
				{
					throw new Exception("This should not have happened. The distributed matrix is not created correctly.");
				}

				if (localRowIndices.TryGetValue(nodeID, out int localRowIdx))
				{
					LocalMatrices[nodeID].SetEntryRespectingPattern(localRowIdx, localColIdx, value); // Do this in all instances.
				}
			}

			// If we reached this line, then the entry (rowIdx, colIdx) is not explicitly stored (structural zero)
			if (value != 0.0)
			{
				throw new SparsityPatternModifiedException(
					$"The entry ({rowIdx}, {colIdx}) is a structural zero and cannot be changed");
			}
		}

		IMatrix IMatrixView.Transpose() => Transpose();

		public DistributedOverlappingMatrix<TMatrix> Transpose()
		{
			var transpose = new DistributedOverlappingMatrix<TMatrix>(Indexer);
			Environment.DoPerNode(nodeID => transpose.LocalMatrices[nodeID] = (TMatrix)this.LocalMatrices[nodeID].Transpose());
			return transpose;
		}

		private void CreateGlobalIndexerIfMissing()
		{
			if (globalIndexer == null)
			{
				lock (globalIndexer)
				{
					if (globalIndexer == null) // in case another thread created it before this thread got the lock
					{
						globalIndexer = new GlobalIndexer(Indexer);
					}
				}
			}
		}
	}
}
