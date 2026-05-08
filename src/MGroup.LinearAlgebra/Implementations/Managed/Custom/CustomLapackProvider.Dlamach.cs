namespace MGroup.LinearAlgebra.Implementations.Managed.Custom
{
	using System;
	using System.Diagnostics;

	public partial class CustomLapackProvider : ILapackProvider
	{
		private const double MachineEpsilon = 2.2204460492503131e-16;
		private const double SafeMinimum = 2.2250738585072014e-308;

		private double Dlamch(DlamchParam param)
		{
			return param switch
			{
				DlamchParam.Epsilon => MachineEpsilon,
				DlamchParam.SafeMinimum => SafeMinimum,
				DlamchParam.Base => 2.0,
				DlamchParam.Precision => MachineEpsilon * 2.0,
				DlamchParam.MantissaDigits => 53.0,
				DlamchParam.Rounding => 1.0,
				DlamchParam.MinExponent => -1022.0,
				DlamchParam.UnderflowThreshold => SafeMinimum,
				DlamchParam.MaxExponent => 1023.0,
				DlamchParam.OverflowThreshold => double.MaxValue,
				_ => throw new ArgumentException("Invalid DLAMCH parameter")
			};
		}
	}
}
