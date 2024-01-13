namespace MGroup.LinearAlgebra.Iterative.PreconditionedConjugateGradient
{
	using System;
	using System.Collections.Generic;
	using System.Text;
	using MGroup.LinearAlgebra.Vectors;

	/** The coefficients of a vector's linear combination.
	 * They are the coefficients only from a linear combination.
	 * The result vector is a linear combination of vectors of 2 Krylov subspaces R and P.
	 * Coefficients of Krylov subspace R is in r list and for Krylov subspace P in p list.
	 * Lists can have fewer (even zero) elements than Krylov subspace vectors. */
	internal class BlockVectorOperator
	{
		private readonly List<double> r;
		private readonly List<double> p;

		public List<double> R { get => r; }
		public List<double> P { get => p; }

		public BlockVectorOperator(int n)
		{
			r = new List<double>(n);
			p = new List<double>(n);
		}

		private static void Xupdate(List<double> x, double a, List<double> p)
		{
			int to2 = p.Count, to1 = Math.Min(to2, x.Count);
			for (var i = 0; i < to1; ++i)
				x[i] += a * p[i];
			for (var i = to1; i < to2; ++i)
				x.Add(a * p[i]);
		}

		/** It is the operation x += a * p.
		 * x and p are coefficients of Krylov subspace vectors linear combination. */
		public void UpdateX(double a, BlockVectorOperator p)
		{
			Xupdate(r, a, p.r);
			Xupdate(this.p, a, p.p);
		}

		private static void Rupdate(List<double> r, double a, List<double> p)
		{
			int to2 = p.Count, to1 = Math.Min(to2, r.Count - 1);
			for (var i = 0; i < to1; ++i)
				r[i + 1] -= a * p[i];
			if (to1 == -1) { r.Add(0); ++to1; }
			for (var i = to1; i < to2; ++i)
				r.Add(-a * p[i]);
		}

		/** It is the operation r -= a * A * p.
		 * x and p are coefficients of Krylov subspace vectors linear combination. */
		public void UpdateR(double a, BlockVectorOperator p)
		{
			Rupdate(r, a, p.r);
			Rupdate(this.p, a, p.p);
		}

		private static void Pupdate(List<double> r, double a, List<double> p)
		{
			for (var i = 0; i < p.Count; ++i) p[i] *= a;
			int to2 = r.Count, to1 = Math.Min(to2, p.Count);
			for (var i = 0; i < to1; ++i)
				p[i] += r[i];
			for (var i = to1; i < to2; ++i)
				p.Add(r[i]);
		}

		/** It is the operation p = r + a * p.
		 * r and p are coefficients of Krylov subspace vectors linear combination. */
		public void UpdateP(BlockVectorOperator r, double a)
		{
			Pupdate(r.r, a, this.r);
			Pupdate(r.p, a, p);
		}


		/** Calculates the linear combination of a vector.
		 * The operation is x = p_coefficints * P + r_coefficients * R
		 * <param name="inf">Block information with Krylov subspaces</param>
		 * <param name="xc">Coefficients of Krylov subspaces vectors</param>
		 * <return>The result vector</return>
		 */
		public IVector EvaluateVector(IVector[] residualKernels, IVector[] directionKernels)
		{
			var x = residualKernels[0].CreateZeroVectorWithSameFormat();
			for (var i = 0; i < r.Count; ++i)
				x.AxpyIntoThis(residualKernels[i], r[i]);
			for (var i = 0; i < p.Count; ++i)
				x.AxpyIntoThis(directionKernels[i], p[i]);
			return x;
		}

		private static double Square1(IList<double> RR, IList<double> r)
		{
			double a = 0;
			for (var i = 0; i < r.Count; ++i)
				a += r[i] * r[i] * RR[2 * i];
			return a;
		}

		private static double Square2(IList<double> RR, IList<double> r)
		{
			double a = 0;
			for (var i = 0; i < r.Count - 1; ++i)
				for (var j = i + 1; j < r.Count; j++)
					a += 2 * r[i] * r[j] * RR[i + j];
			return a;
		}

		private static double Square3(IList<double> RP, IList<double> r, List<double> p)
		{
			double a = 0;
			for (var i = 0; i < r.Count; ++i)
				for (var j = 0; j < p.Count; j++)
					a += 2 * r[i] * p[j] * RP[i + j];
			return a;
		}

		/** Gives the square of vector r.
		 * Actually is not r * r but r * M * r where M is the inverse preconditioner matrix.
		 * Vector r is a linear combination of vectors of Krylov subspace R and Krylov subspace P.
		 * So r * M * r is the sum of a banch of dot products between vectors of Krylov subspace R, of Krylov subspace P and between them.
		 * But we have already the dot products between vectors of Krylov subspaces in vectors RR, PP, and RP.
		 * So we use coefficients of linear combination (described above), stored in vectors rr and rp
		 * and the vectors RR, PP and RP.
		 */
		public double Square(double[] residualSandwiches, double[] directionSandwiches, double[] residualDirectionSandwiches) =>
			Square1(residualSandwiches, r) + Square1(directionSandwiches, p) + Square2(residualSandwiches, r) + Square2(directionSandwiches, p) + 
			Square3(residualDirectionSandwiches, r, p);

		private static double Sandwich1(IList<double> RR, List<double> r)
		{
			double a = 0;
			for (var i = 0; i < r.Count; ++i)
				a += r[i] * r[i] * RR[2 * i + 1];
			return a;
		}

		private static double Sandwich2(IList<double> RR, IList<double> r)
		{
			double a = 0;
			for (var i = 0; i < r.Count - 1; ++i)
				for (var j = i + 1; j < r.Count; j++)
					a += 2 * r[i] * r[j] * RR[i + j + 1];
			return a;
		}

		private static double Sandwich3(IList<double> RP, IList<double> r, IList<double> p)
		{
			double a = 0;
			for (var i = 0; i < r.Count; ++i)
				for (var j = 0; j < p.Count; j++)
					a += 2 * r[i] * p[j] * RP[i + j + 1];
			return a;
		}

		/** Gives the sandwich product of matrix A and vector p.
		 * Actually is not p * A * p but p * M * A * M * p where M is the inverse preconditioner matrix.
		 * Vector p is a linear combination of vectors of Krylov subspace R and Krylov subspace P.
		 * So p * M * A * M * p is the sum of a banch of dot products between vectors of Krylov subspace R, of Krylov subspace P and between them.
		 * But we have already the dot products between vectors of Krylov subspaces in vectors RR, PP, and RP.
		 * So we use coefficients of linear combination (described above), stored in vectors pr and pp
		 * and the vectors RR, PP and RP.
		 */
		public double Sandwich(double[] residualSandwiches, double[] directionSandwiches, double[] residualDirectionSandwiches) =>
			Sandwich1(residualSandwiches, r) + Sandwich1(directionSandwiches, p) + Sandwich2(residualSandwiches, r) + Sandwich2(directionSandwiches, p) + 
			Sandwich3(residualDirectionSandwiches, r, p);
	}

}
