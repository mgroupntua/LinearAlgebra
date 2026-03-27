namespace MGroup.LinearAlgebra.Implementations.Managed
{
	using System;
	using System.Diagnostics;

	public partial class ManagedBlasProvider : IBlasProvider
	{
		public double Dnrm2(int n, double[] x, int offsetX, int incX)
		{
			Debug.Assert(x != null);
			Debug.Assert(offsetX >= 0);
			Debug.Assert(incX != 0);

			if (n < 1 || incX < 1)
			{
				return 0.0;
			}

			if (n == 1)
			{
				return Math.Abs(x[offsetX]);
			}

			double scale = 0.0;
			double sumsq = 1.0;

			int index = offsetX;

			for (int i = 0; i < n; i++)
			{
				double value = x[index];

				if (value != 0.0)
				{
					double absValue = Math.Abs(value);

					if (scale < absValue)
					{
						double ratio = scale / absValue;
						sumsq = 1.0 + sumsq * ratio * ratio;
						scale = absValue;
					}
					else
					{
						double ratio = absValue / scale;
						sumsq += ratio * ratio;
					}
				}

				index += incX;
			}

			return scale * Math.Sqrt(sumsq);
		}
	}
}
