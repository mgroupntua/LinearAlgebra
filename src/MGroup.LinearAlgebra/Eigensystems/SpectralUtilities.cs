namespace MGroup.LinearAlgebra.Eigensystems
{
	using System;
	using System.Collections;
	using System.Collections.Generic;
	using System.Linq;
	using System.Text;

	using MGroup.LinearAlgebra.Matrices;
	using MGroup.LinearAlgebra.Vectors;

	public static class SpectralUtilities
	{
		internal static (Vector valuesSorted, Matrix vectorsSorted) SortSingularValues(Vector values, Matrix vectors, 
			bool descending)
		{
			Vector valuesSorted = values.Copy();
			int[] mapNewToOld = Enumerable.Range(0, values.Length).ToArray();
			if (descending)
			{
				Array.Sort(valuesSorted.RawData, mapNewToOld, new DescendingComparer());
			}
			else
			{
				Array.Sort(valuesSorted.RawData, mapNewToOld);
			}

			Matrix vectorsSorted = vectors.ReorderColumns(mapNewToOld, oldToNew: false);
			return (valuesSorted, vectorsSorted);
		}

		private class DescendingComparer : IComparer
		{
			public int Compare(object x, object y) => Math.Sign((double)y - (double)x); // reverse order
		}
	}
}
