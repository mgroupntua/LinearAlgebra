using System;
using System.Collections.Generic;
using System.Text;

namespace MGroup.LinearAlgebra.Tests.Utilities
{
	public static class ArrayUtilities
	{
		public static bool AreEqual(int[] expected, int[] computed)
		{
			if (expected.Length != computed.Length)
			{
				return false;
			}
			for (int i = 0; i < expected.Length; ++i)
			{
				if (expected[i] != computed[i])
				{
					return false;
				}
			}
			return true;
		}

		public static bool AreEqual(int[,] expected, int[,] computed)
		{
			if ((expected.GetLength(0) != computed.GetLength(0)) || (expected.GetLength(1) != computed.GetLength(1)))
			{
				return false;
			}
			for (int i = 0; i < expected.GetLength(0); ++i)
			{
				for (int j = 0; j < expected.GetLength(1); ++j)
				{
					if (expected[i, j] != computed[i, j])
					{
						return false;
					}
				}
			}
			return true;
		}
	}
}
