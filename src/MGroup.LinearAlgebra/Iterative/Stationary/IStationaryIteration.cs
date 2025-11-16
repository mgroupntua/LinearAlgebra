namespace MGroup.LinearAlgebra.Iterative.Stationary
{
	using System;
	using System.Collections.Generic;
	using System.Text;

	using MGroup.LinearAlgebra.Matrices;
	using MGroup.LinearAlgebra.Vectors;

	public interface IStationaryIteration : ISettingsCopiable<IStationaryIteration>
	{
		public string Name { get; }

		public void Execute(Vector rhs, Vector solution);

		/// <summary>
		/// Use this method, so that this <see cref="IStationaryIteration"/> will try to reuse work done by
		/// <paramref name="other"/> whenever possible.
		/// </summary>
		/// <param name="other"></param>
		public void LinkWith(IStationaryIteration other);

		public void UpdateMatrix(IReadOnlyMatrix matrix, bool isPatternModified);

	}
}
