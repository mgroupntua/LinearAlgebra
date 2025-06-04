using System;
using System.Collections.Generic;
using System.Text;

namespace MGroup.LinearAlgebra.Distributed.Utilities
{
    public static class DataStructureExtensions
    {
		/// <summary>
		/// Creates a new Dictionary with the same keys as the <paramref name="original"/>, but the values are converted
		/// using <paramref name="transform"/>.
		/// </summary>
		public static Dictionary<TKey, TVal2> MapDictionary<TKey, TVal1, TVal2>(
			this Dictionary<TKey, TVal1> original, Func<TVal1, TVal2> transform)
		{
			var result = new Dictionary<TKey, TVal2>();
			foreach (var pair in original)
			{
				result[pair.Key] = transform(pair.Value);
			}

			return result;
		}
    }
}
