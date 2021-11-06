using System;
using System.Collections.Generic;
using System.Text;
using System.Linq;
using MGroup.LinearAlgebra.Vectors;
using Xunit;

namespace MGroup.LinearAlgebra.Tests.Vectors
{
	public static class VectorExtensionsTests
	{
		[Fact]
		private static void TestFind()
		{
			var vector1 = Vector.CreateFromArray(new double[] { 1, 3, 5, 7, 9, 11, 13, 15, 17, 19 });
			Assert.Equal(new int[] { 6 }, vector1.Find(x => x == 13));

			var vector2 = Vector.CreateFromArray(new double[] { 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 });
			Assert.Equal(new int[] { 3 }, vector2.Find(x => Math.Abs(x - 0.3) < 0.001));
		}
	}
}
