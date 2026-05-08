namespace MGroup.LinearAlgebra.Implementations.Managed
{
	using System;
	using System.Diagnostics;

	public partial class ManagedLapackProvider : ILapackProvider
	{
		private static bool LSame(char a, char b)
		{
			return char.ToUpperInvariant(a) == char.ToUpperInvariant(b);
		}

		private static bool LSame(string a, string b)
		{
			Debug.Assert(!string.IsNullOrEmpty(a));
			Debug.Assert(!string.IsNullOrEmpty(b));

			return char.ToUpperInvariant(a[0]) == char.ToUpperInvariant(b[0]);
		}
	}
}
