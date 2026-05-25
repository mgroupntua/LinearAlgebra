namespace MGroup.LinearAlgebra.Implementations.Managed.Custom
{
	using System;
	using System.Diagnostics;

	public partial class CustomLapackProvider : ILapackProvider
	{
		/// <summary>
		/// Hardware compliance test. It checks:
		/// Division by zero → Inf.
		/// Negative division → -Inf.
		/// NaN propagation works correctly.
		/// </summary>
		private int IeeeCheck(int ispec)
		{
			double zero = 0.0;
			double one = 1.0;

			// Check infinity behavior
			double posInf = one / zero;
			if (!double.IsInfinity(posInf) || posInf <= one)
			{
				return 0;
			}

			double negInf = -one / zero;
			if (!double.IsNegativeInfinity(negInf))
			{
				return 0;
			}

			double nan = zero / zero;

			if (ispec == 0)
			{
				return 1;
			}

			// Check NaN propagation
			if (!double.IsNaN(nan))
			{
				return 0;
			}

			if (nan == nan)
			{
				return 0; // NaN must not equal itself
			}

			return 1;
		}
	}
}
