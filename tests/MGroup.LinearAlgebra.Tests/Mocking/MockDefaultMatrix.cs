namespace MGroup.LinearAlgebra.Tests.Mocking
{
	using System;

	using MGroup.LinearAlgebra.Matrices;

	[Serializable]
	internal class MockDefaultMatrix : DefaultMatrix
	{
		private readonly double[,] data;

		public MockDefaultMatrix(int numRows, int numCols)
		{
			this.data = new double[numRows, numCols];
		}

		public MockDefaultMatrix(double[,] matrix)
		{
			this.data = new double[matrix.GetLength(0), matrix.GetLength(1)];
			Array.Copy(matrix, this.data, matrix.Length);
		}

		public override int NumColumns => data.GetLength(1);

		public override int NumRows => data.GetLength(0);

		public override double this[int rowIdx, int colIdx]
		{
			get => data[rowIdx, colIdx];
			set => data[rowIdx, colIdx] = value;
		}

		public override void Clear() => Array.Clear(data);

		public override IMatrix CreateZeroMatrixWithSameFormat() => new MockDefaultMatrix(NumRows, NumColumns);

		public override bool HasSameFormat(IReadOnlyMatrix otherMatrix)
		{
			if (otherMatrix is MockDefaultMatrix casted && casted.NumRows == this.NumRows && casted.NumColumns == this.NumColumns)
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
