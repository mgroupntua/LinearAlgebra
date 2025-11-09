namespace MGroup.LinearAlgebra.Tests.Mocking
{
	using System;

	using MGroup.LinearAlgebra.Matrices;

	[Serializable]
	internal class MockValuesBackedMatrix : ValuesBackedMatrix<MockValuesBackedMatrix>
	{
		private readonly double[] values;

		public MockValuesBackedMatrix(int numRows, int numCols)
		{
			this.values = new double[numRows * numCols];
			NumRows = numRows;
			NumColumns = numCols;
		}

		public MockValuesBackedMatrix(int numRows, int numCols, double[] values)
		{
			if (values.Length != numRows * numCols)
			{
				throw new ArgumentException("Invalid dimensions");
			}

			this.values = new double[values.Length];
			Array.Copy(values, this.values, values.Length);
			NumRows = numRows;
			NumColumns = numCols;

		}

		public override int NumColumns { get; }

		public override int NumRows { get; }

		public override double[] RawValues => values;

		public override double this[int rowIdx, int colIdx]
		{
			get => values[rowIdx * NumColumns + colIdx];
			set => values[rowIdx * NumColumns + colIdx] = value;
		}

		public override MockValuesBackedMatrix CopyAsSameType(bool copyIndexingData)
			=> new MockValuesBackedMatrix(NumRows, NumColumns, values);
		
		public override MockValuesBackedMatrix CreateZeroMatrixSame()
			=> new MockValuesBackedMatrix(NumRows, NumColumns);

		public override bool HasSameFormat(MockValuesBackedMatrix other)
		{
			if (other.NumRows == this.NumRows && other.NumColumns == this.NumColumns)
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
