namespace MGroup.LinearAlgebra.AlgebraicMultiGrid
{
	using System;
	using System.Collections.Generic;

	using MGroup.LinearAlgebra.Iterative;
	using MGroup.LinearAlgebra.Iterative.Stationary;
	using MGroup.LinearAlgebra.Matrices;
	using MGroup.LinearAlgebra.Vectors;

	public class MultigridLevelSmoothing : ISettingsCopiable<MultigridLevelSmoothing>
	{
		private readonly List<(IStationaryIteration stationaryIteration, int numApplications)> preSmoothers;
		private readonly List<(IStationaryIteration stationaryIteration, int numApplications)> postSmoothers;

		public MultigridLevelSmoothing()
		{
			preSmoothers = new List<(IStationaryIteration, int)>();
			postSmoothers = new List<(IStationaryIteration, int)>();
		}

		public MultigridLevelSmoothing AddPreSmoother(IStationaryIteration stationaryIteration, int numApplications)
		{
			LinkWithExistingSmoothers(stationaryIteration);
			preSmoothers.Add((stationaryIteration, numApplications));
			return this;
		}

		public MultigridLevelSmoothing AddPostSmoother(IStationaryIteration stationaryIteration, int numApplications)
		{
			LinkWithExistingSmoothers(stationaryIteration);
			postSmoothers.Add((stationaryIteration, numApplications));
			return this;
		}

		public void ApplyPreSmoothers(IReadOnlyVector rhs, IVector solution) => ApplySmoothers(preSmoothers, rhs, solution);

		public void ApplyPostSmoothers(IReadOnlyVector rhs, IVector solution) => ApplySmoothers(postSmoothers, rhs, solution);

		public MultigridLevelSmoothing CopyWithInitialSettings()
		{
			var clone = new MultigridLevelSmoothing();
			foreach ((IStationaryIteration stationaryIteration, int numApplications) in preSmoothers)
			{
				clone.AddPreSmoother(stationaryIteration.CopyWithInitialSettings(), numApplications);
			}

			foreach ((IStationaryIteration stationaryIteration, int numApplications) in postSmoothers)
			{
				clone.AddPostSmoother(stationaryIteration.CopyWithInitialSettings(), numApplications);
			}

			return clone;
		}

		public MultigridLevelSmoothing SetPostSmoothersSameAsPreSmoothers()
		{
			postSmoothers.Clear();
			postSmoothers.AddRange(preSmoothers);
			return this;
		}

		public void UpdateMatrix(IReadOnlyMatrix matrix, bool isPatternModified)
		{
			CheckSmoothers();
			foreach ((IStationaryIteration stationaryIteration, _) in preSmoothers)
			{
				stationaryIteration.UpdateMatrix(matrix, isPatternModified);
			}

			foreach ((IStationaryIteration stationaryIteration, _) in postSmoothers)
			{
				stationaryIteration.UpdateMatrix(matrix, isPatternModified);
			}
		}

		private static void ApplySmoothers(List<(IStationaryIteration stationaryIteration, int numApplications)> smoothers,
			IReadOnlyVector rhs, IVector solution)
		{
			var rhsDense = (Vector)rhs;
			var solutionDense = (Vector)solution;
			foreach ((IStationaryIteration stationaryIteration, int numApplications) in smoothers)
			{
				for (var t = 0; t < numApplications; t++)
				{
					stationaryIteration.Execute(rhsDense, solutionDense);
				}
			}
		}

		private void CheckSmoothers()
		{
			if (preSmoothers.Count == 0)
			{
				throw new InvalidOperationException("There are no pre-smoothers defined");
			}

			if (postSmoothers.Count == 0)
			{
				throw new InvalidOperationException("There are no post-smoothers defined");
			}
		}

		private void LinkWithExistingSmoothers(IStationaryIteration newIteration)
		{
			foreach ((IStationaryIteration existingIteration, _) in preSmoothers)
			{
				if (newIteration != existingIteration)
				{
					newIteration.LinkWith(existingIteration);
				}
			}

			foreach ((IStationaryIteration existingIteration, _) in postSmoothers)
			{
				if (newIteration != existingIteration)
				{
					newIteration.LinkWith(existingIteration);
				}
			}
		}
	}
}
