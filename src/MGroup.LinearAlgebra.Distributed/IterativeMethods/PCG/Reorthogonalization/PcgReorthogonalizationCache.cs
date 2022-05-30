using System;
using System.Collections.Generic;
using System.Diagnostics;
using MGroup.LinearAlgebra.Vectors;
using MGroup.MSolve.Solution.LinearSystem;

namespace MGroup.LinearAlgebra.Distributed.IterativeMethods.PCG.Reorthogonalization
{
	/// <summary>
	/// Manages the insertion and removal of PCG direction vectors and related data, that will be used for reorthogonalization.
	/// Authors: Serafeim Bakalakos
	/// </summary>
	public class PcgReorthogonalizationCache
	{
		/// <summary>
		/// The conjugate direction vectors stored so far.
		/// </summary>
		public List<IGlobalVector> Directions { get; } = new List<IGlobalVector>();

		/// <summary>
		/// The products <see cref="Directions"/> * systemMatrix * <see cref="Directions"/> stored so far.
		/// </summary>
		public List<double> DirectionsTimesMatrixTimesDirections { get; } = new List<double>();

		/// <summary>
		/// The index into <see cref="Directions"/> of the first vector of each generation. A generation is defined as all
		/// direction vectors (and related data) stored when solving the linear system for a specific rhs. Thus solving for 2
		/// consecutive rhs vectors will generate direction vectors belonging to 2 generations.
		/// </summary>
		public List<int> GenerationStartIndices { get; } = new List<int>();

		/// <summary>
		/// The products systemMatrix * <see cref="Directions"/> stored so far.
		/// </summary>
		public List<IGlobalVector> MatrixTimesDirections { get; } = new List<IGlobalVector>();

		public bool AreAllDirectionsConjugate(double tolerance)
		{
			int numVectors = Directions.Count;

			// Examine the newest direction first, since it is the most probable to be affected by error build-up
			for (int i = numVectors - 1; i >= 1; --i) 
			{
				for (int j = i - 1; j >= 0; --j)
				{
					// Conjugate if di * A * dj = 0.
					double dot = Directions[i].DotProduct(MatrixTimesDirections[j]);
					if (dot > tolerance)
					{
						return false;
					}
				}
			}

			return true;
		}

		public List<double> CalcMaxConjugacyFactorPerIteration()
		{
			int numVectors = Directions.Count;
			var factors = new List<double>(numVectors);
			factors.Add(double.NaN); // No other vector yet
			for (int i = 1; i < numVectors; ++i)
			{
				double max = -1.0;
				for (int j = i - 1; j >= 0; --j)
				{
					double dot = Math.Abs(Directions[i].DotProduct(MatrixTimesDirections[j]));
					if (dot > max)
					{
						max = dot;
					}
				}
				factors.Add(max);
			}
			return factors;
		}

		public void Clear()
		{
			Directions.Clear();
			MatrixTimesDirections.Clear();
			DirectionsTimesMatrixTimesDirections.Clear();
		}

		/// <summary>
		/// Discards the direction vectors and any corresponding data of the newest PCG iterations.
		/// </summary>
		/// <param name="numOldVectorsToRemove">
		/// The number of the newest entries (direction vectors and corresponding data) to discard. If it exceeds the number of
		/// entries currently stored, they will all be discarded without throwing any exceptions.
		/// </param>
		public void RemoveNewDirectionVectorData(int numNewVectorsToRemove)
		{
			if (numNewVectorsToRemove > Directions.Count)
			{
				Directions.Clear();
				MatrixTimesDirections.Clear();
				DirectionsTimesMatrixTimesDirections.Clear();
			}
			else
			{
				int start = Directions.Count - numNewVectorsToRemove;
				Directions.RemoveRange(start, numNewVectorsToRemove);
				MatrixTimesDirections.RemoveRange(start, numNewVectorsToRemove);
				DirectionsTimesMatrixTimesDirections.RemoveRange(start, numNewVectorsToRemove);
			}
		}

		/// <summary>
		/// Discards the direction vectors and any corresponding data of the oldest PCG iterations.
		/// </summary>
		/// <param name="numOldVectorsToRemove">
		/// The number of the oldest entries (direction vectors and corresponding data) to discard. If it exceeds the number of
		/// entries currently stored, they will all be discarded without throwing any exceptions.
		/// </param>
		public void RemoveOldDirectionVectorData(int numOldVectorsToRemove)
		{
			if (numOldVectorsToRemove > Directions.Count)
			{
				Directions.Clear();
				MatrixTimesDirections.Clear();
				DirectionsTimesMatrixTimesDirections.Clear();
			}
			else
			{
				Directions.RemoveRange(0, numOldVectorsToRemove);
				MatrixTimesDirections.RemoveRange(0, numOldVectorsToRemove);
				DirectionsTimesMatrixTimesDirections.RemoveRange(0, numOldVectorsToRemove);
			}
		}

		public void StartGeneration()
		{
			GenerationStartIndices.Add(Directions.Count);
		}

		/// <summary>
		/// Stores a new direction vector and other related data. The new entries will be regarded as latest.
		/// </summary>
		/// <param name="pcg">The Preconditioned Conjugate Gradient Aglorithm that uses this object.</param>
		public void StoreDirectionData(ReorthogonalizedPcg pcg)
		{
			Directions.Add(pcg.Direction.Copy());
			MatrixTimesDirections.Add(pcg.MatrixTimesDirection.Copy());
			DirectionsTimesMatrixTimesDirections.Add(pcg.DirectionTimesMatrixTimesDirection);
		}
	}
}
