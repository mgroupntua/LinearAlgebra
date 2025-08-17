namespace MGroup.LinearAlgebra.Commons
{
	using System;
	using System.Diagnostics;

	using MGroup.LinearAlgebra.Exceptions;

	public static class PerformanceWarnings
	{
		[Conditional("RELEASE")]
		public static void ProhibitPerformanceBottlenecks()
		{
			if (LibrarySettings.ThrowExceptionOnKnownPerformanceBottlenecksInReleaseBuilds)
			{
				throw new PerformanceBottleneckException(
					"Potential performance bottleneck due to accessing all entries of a potentially sparse matrix or vector.");
			}
		}

		[Conditional("DEBUG")]
		public static void WarnAboutPerformanceBottlenecks()
		{
			Debug.WriteLine(
				"Potential performance bottleneck due to accessing all entries of a potentially sparse matrix or vector, at :"
				+ Environment.StackTrace);
		}
	}
}
