namespace MGroup.LinearAlgebra.Iterative.PreconditionedConjugateGradient
{
	using System;
	using System.Collections.Generic;
	using System.Text;
	using MGroup.LinearAlgebra.Vectors;

	/// <summary>
	/// This class handles linear combination calculation pertaining to the block operations of the Krylov subspaces R and P.
	/// </summary>
	internal class BlockVectorOperator
	{
		private readonly List<double> r;
		private readonly List<double> p;

		/// <summary>
		/// This list contains the linear combination coefficients pertaining to subspace R.
		/// </summary>
		public List<double> R { get => r; }

		/// <summary>
		/// This list contains the linear combination coefficients pertaining to subspace P.
		/// </summary>
		public List<double> P { get => p; }

		/// <summary>
		/// This class handles linear combination calculation pertaining to the block operations of the Krylov subspaces R and P.
		/// </summary>
		/// <param name="n">The dimension of the Krylov subspace which coincides with the number of steps of the block CG algorithm</param>
		public BlockVectorOperator(int n)
		{
			r = new List<double>(n);
			p = new List<double>(n);
		}

		private static void Xupdate(List<double> x, double a, List<double> p)
		{
			int to2 = p.Count, to1 = Math.Min(to2, x.Count);
			for (var i = 0; i < to1; ++i)
			{
				x[i] += a * p[i];
			}
			
			for (var i = to1; i < to2; ++i)
			{
				x.Add(a * p[i]);
			}
		}

		/// <summary>
		/// Calculates x += a * p with x and p being the coefficients of Krylov subspace vectors linear combination./>.
		/// </summary>
		/// <param name="a">The value of scalar a of the operation x += a * p.</param>
		/// <param name="p">The p coefficients contained in a <see cref="BlockVectorOperator"/> object.</param>
		public void UpdateX(double a, BlockVectorOperator p)
		{
			Xupdate(r, a, p.r);
			Xupdate(this.p, a, p.p);
		}

		private static void Rupdate(List<double> r, double a, List<double> p)
		{
			int to2 = p.Count, to1 = Math.Min(to2, r.Count - 1);
			for (var i = 0; i < to1; ++i)
			{
				r[i + 1] -= a * p[i];
			}
			
			if (to1 == -1) 
			{ 
				r.Add(0); 
				++to1; 
			}
			
			for (var i = to1; i < to2; ++i)
			{
				r.Add(-a * p[i]);
			}
		}

		/// <summary>
		/// Calculates r -= a * A * p with a and p being the coefficients of Krylov subspace vectors linear combination./>.
		/// </summary>
		/// <param name="a">The value of scalar a of the operation r -= a * A * p.</param>
		/// <param name="p">The r coefficients contained in a <see cref="BlockVectorOperator"/> object.</param>
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
			{
				p[i] += r[i];
			}
			
			for (var i = to1; i < to2; ++i)
			{
				p.Add(r[i]);
			}
		}

		/// <summary>
		/// Calculates p = r + a * p with r and a being the coefficients of Krylov subspace vectors linear combination./>.
		/// </summary>
		/// <param name="r">The r coefficients contained in a <see cref="BlockVectorOperator"/> object.</param>
		/// <param name="a">The value of scalar a of the operation p = r + a * p.</param>
		public void UpdateP(BlockVectorOperator r, double a)
		{
			Pupdate(r.r, a, this.r);
			Pupdate(r.p, a, p);
		}

		/// <summary>
		/// Calculates x = p_coefficients * P + r_coefficients * R and returns the result of x/>.
		/// </summary>
		/// <param name="residualKernels">The r_coefficients.</param>
		/// <param name="directionKernels">The p_coefficients.</param>
		/// <returns>The result of x = p_coefficients * P + r_coefficients * R</returns>
		public IVector EvaluateVector(IVector[] residualKernels, IVector[] directionKernels)
		{
			var x = residualKernels[0].CreateZeroVectorWithSameFormat();
			for (var i = 0; i < r.Count; ++i)
			{
				x.AxpyIntoThis(residualKernels[i], r[i]);
			}
			
			for (var i = 0; i < p.Count; ++i)
			{
				x.AxpyIntoThis(directionKernels[i], p[i]);
			}
			
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
			{
				for (var j = i + 1; j < r.Count; j++)
				{
					a += 2 * r[i] * r[j] * RR[i + j];
				}
			}
			
			return a;
		}

		private static double Square3(IList<double> RP, IList<double> r, List<double> p)
		{
			double a = 0;
			for (var i = 0; i < r.Count; ++i)
			{
				for (var j = 0; j < p.Count; j++)
				{
					a += 2 * r[i] * p[j] * RP[i + j];
				}
			}
			return a;
		}

		/// <summary>
		/// Calculates r * M * r where M is the inverse preconditioner matrix. Vector r is a linear combination of vectors of Krylov subspace R and Krylov subspace P./>.
		/// </summary>
		/// <param name="residualSandwiches">The RR.</param>
		/// <param name="directionSandwiches">The PP.</param>
		/// <param name="residualDirectionSandwiches">The RP.</param>
		/// <returns>The result of r * M * r</returns>
		/// <remarks>
		/// Vector r is a linear combination of vectors of Krylov subspace R and Krylov subspace P.
		/// So r * M * r is the sum of a banch of dot products between vectors of Krylov subspace R, of Krylov subspace P and between them.
		/// But we have already the dot products between vectors of Krylov subspaces in vectors RR, PP, and RP.
		/// So we use coefficients of linear combination (described above), stored in vectors rr and rp
		/// and the vectors RR, PP and RP.
		/// </remarks>
		public double Square(double[] residualSandwiches, double[] directionSandwiches, double[] residualDirectionSandwiches) =>
			Square1(residualSandwiches, r) + Square1(directionSandwiches, p) + Square2(residualSandwiches, r) + Square2(directionSandwiches, p) + 
			Square3(residualDirectionSandwiches, r, p);

		private static double Sandwich1(IList<double> RR, List<double> r)
		{
			double a = 0;
			for (var i = 0; i < r.Count; ++i)
			{
				a += r[i] * r[i] * RR[2 * i + 1];
			}
			
			return a;
		}

		private static double Sandwich2(IList<double> RR, IList<double> r)
		{
			double a = 0;
			for (var i = 0; i < r.Count - 1; ++i)
			{
				for (var j = i + 1; j < r.Count; j++)
				{
					a += 2 * r[i] * r[j] * RR[i + j + 1];
				}
			}
			
			return a;
		}

		private static double Sandwich3(IList<double> RP, IList<double> r, IList<double> p)
		{
			double a = 0;
			for (var i = 0; i < r.Count; ++i)
			{
				for (var j = 0; j < p.Count; j++)
				{
					a += 2 * r[i] * p[j] * RP[i + j + 1];
				}
			}
			
			return a;
		}

		/// <summary>
		/// Calculates p * M * A * M * p (sandwich product) where M is the inverse preconditioner matrix. Vector p is a linear combination of vectors of Krylov subspace R and Krylov subspace P./>.
		/// </summary>
		/// <param name="residualSandwiches">The RR.</param>
		/// <param name="directionSandwiches">The PP.</param>
		/// <param name="residualDirectionSandwiches">The RP.</param>
		/// <returns>The result of r * M * r</returns>
		/// <remarks>
		/// Vector p is a linear combination of vectors of Krylov subspace R and Krylov subspace P.
		/// So p * M * A * M * p is the sum of a banch of dot products between vectors of Krylov subspace R, of Krylov subspace P and between them.
		/// But we have already the dot products between vectors of Krylov subspaces in vectors RR, PP, and RP.
		/// So we use coefficients of linear combination (described above), stored in vectors pr and pp
		/// and the vectors RR, PP and RP.
		/// </remarks>		
		public double Sandwich(double[] residualSandwiches, double[] directionSandwiches, double[] residualDirectionSandwiches) =>
			Sandwich1(residualSandwiches, r) + Sandwich1(directionSandwiches, p) + Sandwich2(residualSandwiches, r) + Sandwich2(directionSandwiches, p) + 
			Sandwich3(residualDirectionSandwiches, r, p);
	}
}
