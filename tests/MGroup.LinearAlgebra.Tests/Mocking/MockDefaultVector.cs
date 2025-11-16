namespace MGroup.LinearAlgebra.Tests.Mocking
{
	using System;

	using MGroup.LinearAlgebra.Vectors;

	internal class MockDefaultVector : DefaultVector
	{
		private readonly double[] data;

		public MockDefaultVector(int length)
		{
			this.data = new double[length];
		}

		public MockDefaultVector(double[] values)
		{
			this.data = new double[values.Length];
			Array.Copy(values, this.data, values.Length);
		}

		public override int Length => data.Length;

		public override double this[int index]
		{
			get => data[index];
			set => data[index] = value;
		}

		public override void Clear() => Array.Clear(data);

		public override IVector CreateZeroVectorWithSameFormat() => new MockDefaultVector(Length);

		public override bool HasSameFormat(IReadOnlyVector other)
		{
			if (other is MockDefaultVector casted && casted.Length == this.Length)
			{
				return true;
			}
			else
			{
				return false;
			}
		}
	}
}
