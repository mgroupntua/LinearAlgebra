namespace MGroup.LinearAlgebra.Tests.Utilities
{
	using MGroup.LinearAlgebra.Commons;
	using MGroup.LinearAlgebra.Matrices;
	using MGroup.LinearAlgebra.SchurComplements.IntegerMatrices;
	using MGroup.LinearAlgebra.Vectors;

	using Xunit;

	/// <summary>
	/// Compares scalars, vectors and matrices.
	/// </summary>
	public class MatrixComparer
	{
		private readonly ValueComparer valueComparer;

		public MatrixComparer(double tolerance = 1e-13)
		{
			this.valueComparer = new ValueComparer(tolerance);
		}

		public bool AreEqual(double a, double b) => valueComparer.AreEqual(a, b);

		public bool AreEqual(int[] a, int[] b)
		{
			int n = a.Length;
			if (b.Length != n)
			{
				return false;
			}

			for (int i = 0; i < n; ++i)
			{
				if (a[i] != b[i])
				{
					return false;
				}
			}

			return true;
		}

		public bool AreEqual(double[] a, double[] b)
		{
			int n = a.Length;
			if (b.Length != n)
			{
				return false;
			}

			for (int i = 0; i < n; ++i)
			{
				if (!valueComparer.AreEqual(a[i], b[i]))
				{
					return false;
				}
			}

			return true;
		}

		public bool AreEqual(double[,] a, double[,] b)
		{
			int m = a.GetLength(0);
			int n = a.GetLength(1);
			if ((b.GetLength(0) != m) || (b.GetLength(1) != n))
			{
				return false;
			}

			for (int i = 0; i < m; ++i)
			{
				for (int j = 0; j < n; ++j)
				{
					if (!valueComparer.AreEqual(a[i, j], b[i, j]))
					{
						return false;
					}
				}
			}

			return true;
		}

		public bool AreEqual(IIndexable1D a, IIndexable1D b)
		{
			int n = a.Length;
			if (b.Length != n)
			{
				return false;
			}

			for (int i = 0; i < n; ++i)
			{
				if (!valueComparer.AreEqual(a[i], b[i]))
				{
					return false;
				}
			}

			return true;
		}

		public bool AreEqual(double[] a, IIndexable1D b)
		{
			int n = a.Length;
			if (b.Length != n)
			{
				return false;
			}

			for (int i = 0; i < n; ++i)
			{
				if (!valueComparer.AreEqual(a[i], b[i]))
				{
					return false;
				}
			}

			return true;
		}

		public bool AreEqual(IIndexable1D a, double[] b) => AreEqual(b, a);

		public bool AreEqual(IIndexable2D a, IIndexable2D b)
		{
			int m = a.NumRows;
			int n = a.NumColumns;
			if ((b.NumRows != m) || (b.NumColumns != n))
			{
				return false;
			}

			for (int i = 0; i < m; ++i)
			{
				for (int j = 0; j < n; ++j)
				{
					if (!valueComparer.AreEqual(a[i, j], b[i, j]))
					{
						return false;
					}
				}
			}

			return true;
		}

		public bool AreEqual(double[,] a, IIndexable2D b)
		{
			int m = a.GetLength(0);
			int n = a.GetLength(1);
			if ((b.NumRows != m) || (b.NumColumns != n))
			{
				return false;
			}

			for (int i = 0; i < m; ++i)
			{
				for (int j = 0; j < n; ++j)
				{
					if (!valueComparer.AreEqual(a[i, j], b[i, j]))
					{
						return false;
					}
				}
			}

			return true;
		}

		public bool AreEqual(IIndexable2D a, double[,] b) => AreEqual(b, a);


		public bool AreEqual(int[,] a, IIndexableInt2D b)
		{
			int m = a.GetLength(0);
			int n = a.GetLength(1);
			if ((b.NumRows != m) || (b.NumColumns != n))
			{
				return false;
			}

			for (int i = 0; i < m; ++i)
			{
				for (int j = 0; j < n; ++j)
				{
					if (a[i, j] != b[i, j])
					{
						return false;
					}
				}
			}

			return true;
		}

		public bool AreEqual(IIndexableInt2D a, int[,] b) => AreEqual(b, a);

		public void AssertEqual(double a, double b) => Assert.True(AreEqual(a, b), $"a={a}, b={b}");

		public void AssertEqual(int[] a, int[] b) => Assert.True(AreEqual(a, b));

		public void AssertEqual(double[] a, double[] b) => Assert.True(AreEqual(a, b));

		public void AssertEqual(double[,] a, double[,] b) => Assert.True(AreEqual(a, b));

		public void AssertEqual(IIndexable1D a, IIndexable1D b) => Assert.True(AreEqual(a, b));

		public void AssertEqual(double[] a, IIndexable1D b) => Assert.True(AreEqual(a, b));

		public void AssertEqual(IIndexable1D a, double[] b) => Assert.True(AreEqual(a, b));

		public void AssertEqual(IIndexable2D a, IIndexable2D b) => Assert.True(AreEqual(a, b));

		public void AssertEqual(double[,] a, IIndexable2D b) => Assert.True(AreEqual(a, b));

		public void AssertEqual(IIndexable2D a, double[,] b) => Assert.True(AreEqual(a, b));

		public void AssertEqual(int[,] a, IIndexableInt2D b) => Assert.True(AreEqual(a, b));

		public void AssertEqual(IIndexableInt2D a, int[,] b) => Assert.True(AreEqual(a, b));
	}
}
