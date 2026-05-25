namespace MGroup.LinearAlgebra.Implementations.Managed.Custom
{
	using System;
	using System.Diagnostics;

	public partial class CustomLapackProvider : ILapackProvider
	{
		private int Iparmq(IlaenvSpec spec, int n, int ilo, int ihi, int lwork)
		{
			int nh = ihi - ilo + 1;

			int ns = 2;
			if (nh >= 30) ns = 4;
			if (nh >= 60) ns = 10;
			if (nh >= 150) ns = Math.Max(10, nh / (int)Math.Log(nh, 2));
			if (nh >= 590) ns = 64;
			if (nh >= 3000) ns = 128;
			if (nh >= 6000) ns = 256;

			ns = Math.Max(2, ns - (ns % 2)); // make even

			return spec switch
			{
				IlaenvSpec.IpArmq => 75,              // ISPEC 12
				IlaenvSpec.MinimumBlockSize => 14,    // ISPEC 14
				IlaenvSpec.NumberOfShifts => ns,      // ISPEC 15
				IlaenvSpec.MinimumColumnDimension => nh > 500 ? (3 * ns / 2) : ns, // ISPEC 13
				IlaenvSpec.CrossoverPoint => ns >= 14 ? 2 : 0, // ISPEC 16
				_ => -1
			};
		}
	}
}
