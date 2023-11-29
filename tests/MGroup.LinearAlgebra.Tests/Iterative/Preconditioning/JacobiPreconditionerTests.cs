namespace MGroup.LinearAlgebra.Tests.Iterative.Preconditioning
{
	using System;
	using System.Collections.Generic;
	using System.Linq;
	using System.Text;
	using System.Threading.Tasks;

	using MGroup.LinearAlgebra.Iterative.Preconditioning;
	using MGroup.LinearAlgebra.Matrices;
	using MGroup.LinearAlgebra.Tests.TestData;
	using MGroup.LinearAlgebra.Tests.Utilities;
	using MGroup.LinearAlgebra.Vectors;

	using Xunit;

	public class JacobiPreconditionerTests
	{
		private static readonly MatrixComparer comparer = new MatrixComparer(1E-10);


		[Theory]
		[InlineData(true)]
		[InlineData(false)]
		public static void TestSingleIteration(bool preinvert)
		{
			var matrix = Matrix.CreateFromArray(SquareSingular10by10.Matrix);
			var factory = new JacobiPreconditioner.Factory();
			factory.PreInvert = preinvert;
			var preconditioner = factory.CreatePreconditionerFor(matrix);

			var b = Vector.CreateWithValue(10, 1.0);
			var xExpected = Vector.CreateFromArray(new double[] 
			{ 
				1, 3.000300030003000, 0.25, 0.2, 4, 0.111111111111111, 0.125, 0.142857142857143, 1, 0.5
			});
			var xComputed = Vector.CreateZero(10);

			preconditioner.SolveLinearSystem(b, xComputed);
			comparer.AssertEqual(xExpected, xComputed);
		}
	}
}
