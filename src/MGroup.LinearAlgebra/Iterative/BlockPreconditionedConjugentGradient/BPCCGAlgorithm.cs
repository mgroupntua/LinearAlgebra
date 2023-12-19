namespace MGroup.LinearAlgebra.Iterative.BlockPreconditionedConjugateGradient
{
	using System;
	using System.Collections.Generic;
	using System.Linq;

	using MGroup.LinearAlgebra.Commons;
	using MGroup.LinearAlgebra.Iterative.Preconditioning;
	using MGroup.LinearAlgebra.Iterative.Termination;
	using MGroup.LinearAlgebra.Matrices;
	using MGroup.LinearAlgebra.Vectors;

	public class BPCCGAlgorithm
	{
		/** Creates Krylov subspace from a matrix and a vector.
		 * Normally this function must be run in parallel in multiple kernels.
		 * <param name="linearTransformation">A linear transformation equivaled to matrix-vector product.
		 * It corresponds to the matrix of matrix-vector pair for Krylov subspace.</param>.
		 * <param name="preconditioner">Preconditioner</param>
		 * <param name="vector">The vector of pair matrix - vector for Krylov subsubspace.</param>
		 * <param name="num_of_powers">Number of vectors in Krylov subspace.</param>
		 * <return>The krylov subspace. A list of <paramref name="num_of_powers"/> vectors.
		 * First vector is vector <paramref name="vector"/> and last vector is
		 * (<paramref name="linearTransformation"/> * <paramref name="preconditioner"/>)^<paramref name="num_of_powers"/> * <paramref name="vector"/>.</return>.
		 */
		public static List<IVector> CreateKrylovSubspace(ILinearTransformation linearTransformation, IPreconditioner preconditioner, IVector vector, int num_of_powers)
		{
			List<IVector> krylov = new List<IVector>(num_of_powers);
			krylov.Add(vector);
			for (int i = 1; i <= num_of_powers; ++i)
			{
				IVector v1 = vector.CreateZeroVectorWithSameFormat();
				preconditioner.SolveLinearSystem(krylov[i - 1], v1);
				IVector v2 = vector.CreateZeroVectorWithSameFormat();
				linearTransformation.Multiply(v1, v2);
				krylov.Add(v2);
			}
			return krylov;
		}

		/** Creates Krylov subspace from a matrix, a inverse preconditioner and a vector.
		 * Normally this function must be run in parallel in multiple kernels.
		 * <param name="matrix">The matrix which corresponds to the matrix of matrix-vector pair for Krylov subspace.</param>.
		 * <param name="preconditioner">Preconditioner</param>
		 * <param name="vector">The vector of pair matrix - vector for Krylov subsubspace.</param>
		 * <param name="num_of_powers">Number of vectors in Krylov subspace.</param>
		 * <return>The krylov subspace. A list of n vectors. First vector is vector <paramref name="vector"/>
		 * (<paramref name="matrix"/> * <paramref name="preconditioner"/>)^<paramref name="num_of_powers"/> * <paramref name="vector"/>.</return>.
		 */
		public static List<IVector> CreateKrylovSubspace(IMatrixView matrix, IPreconditioner preconditioner, IVector vector, int num_of_powers)
		{
			List<IVector> krylov = new List<IVector>(num_of_powers);
			krylov.Add(vector);
			for (int i = 1; i <= num_of_powers; ++i)
			{
				IVector v = vector.CreateZeroVectorWithSameFormat();
				preconditioner.SolveLinearSystem(krylov[i - 1], v);
				v = matrix.Multiply(v);
				krylov.Add(v);
			}
			return krylov;
		}

		/** Creates preconditioned dot products with the vectors or 2 Krylov subspaces.
		 * Normally this function must be run in parallel in multiple kernels.
		 * If R(i) is the a vector of R Krylov subspace and P(i) a vector of P Krylov subspace,
		 * this function produces preconditioned dot products R(i) * M * P(i) where M is the inverse preconditioner matrix.
		 * <param name="preconditioner">The preconditioner.</param>
		 * <param name="krylov1">The krylov subspace R. A list of n vectors. First vector is vector (A * M)^0 * r
		 * and last vector is (A * M)^(n-1) * r, where A is the matrix and M is the inverse preconditioner.</param>.
		 * <param name="krylov2">The krylov subspace P. A list of n vectors. First vector is vector (A * M)^0 * p
		 * and last vector is (A * M)^(n-1) * p, where A is the matrix and M is the inverse preconditioner.</param>.
		 * <return>The result vector with the preconditioned dot products. This is the output of this function.</return>
		 */
		public static List<double> CreateKrylovSubspaceDotProducts(IPreconditioner preconditioner, List<IVector> krylov1, List<IVector> krylov2)
		{
			List<double> dots = new List<double>(krylov1.Count + krylov2.Count - 1);
			IVector v = krylov1[0].CreateZeroVectorWithSameFormat();
			preconditioner.SolveLinearSystem(krylov1[0], v);
			for (int i = 0; i < krylov2.Count; ++i)
				dots.Add(v.DotProduct(krylov2[i]));
			preconditioner.SolveLinearSystem(krylov2.Last(), v);
			for (int i = 1; i < krylov1.Count; ++i)
				dots.Add(v.DotProduct(krylov1[i]));
			return dots;
		}


		/** Data calculated, normally, in parallel processing part of code. */
		private struct BlockInfo
		{
			/** Krylov subspace of (A * M)^n * r.
			 * A is the matrix, M is inverse preconditioner matrix and r is the CG residual vector. */
			public List<IVector> R;
			/** Krylov subspace of (A * M)^n * p.
			 * A is the matrix, M is inverse preconditioner matrix and p is the CG conjugate direction vector. */
			public List<IVector> P;
			/** 2n-1 sandwich products of n-vector Krylov subspace (A*M, r).
			 * Sandwich product is r_i * M * r_j. */
			public List<double> RR;
			/** 2n-1 sandwich products of n-vector Krylov subspace (A*M, p).
			 * Sandwich product is p_i * M * p_j. */
			public List<double> PP;
			/** 2n-1 sandwich products of (n-1) vector Krylov subspace (A*M, r) with n-vector Krylov subspace (A*M, p).
			 * Sandwich product is r_i * M * p_j. */
			public List<double> RP;
		};


		/** Create all the information for the first block of algorithm.
 		 * Normally this function must be run in parallel in multiple kernels.
 		 * <param name="linearTransformation">A linear transformation (A * M) * vector where A is the matrix and M is the inverse preconditioner matrix.</param>
 		 * <param name="preconditioner">The preconditioner.</param>
 		 * <param name="r">Vector r of CG.</param>
		 * <param name="n">Number of vectors in Krylov subspace.</param>
 		 * <return>The results of parallel processing.</return>
 		 */
		private static BlockInfo CreateBlockInfo(ILinearTransformation linearTransformation, IPreconditioner preconditioner, IVector r, int n)
		{
			BlockInfo inf = new BlockInfo();
			inf.R = CreateKrylovSubspace(linearTransformation, preconditioner, r, n);
			inf.P = inf.R;
			inf.RR = CreateKrylovSubspaceDotProducts(preconditioner, inf.R, inf.R);
			inf.PP = inf.RR;
			inf.RP = inf.RR;
			return inf;
		}

		/** Create all the information for all the blocks of algorithm except first.
 		 * Normally this function must be run in parallel in multiple kernels.
 		 * <param name="linearTransformation">A linear transformation (A * M) * vector where A is the matrix and M is the inverse preconditioner matrix.</param>
 		 * <param name="preconditioner">The preconditioner.</param>
 		 * <param name="r">Vector r of CG.</param>
 		 * <param name="p">Vector p of CG.</param>
 		 * <param name="n">Number of vectors in Krylov subspace.</param>
 		 * <return>The results of parallel processing.</return>
		 */
		private static BlockInfo CreateBlockInfo(ILinearTransformation linearTransformation, IPreconditioner preconditioner, IVector r, IVector p, int n)
		{
			BlockInfo inf = new BlockInfo();
			inf.R = CreateKrylovSubspace(linearTransformation, preconditioner, r, n);
			inf.P = CreateKrylovSubspace(linearTransformation, preconditioner, p, n + 1);
			inf.RR = CreateKrylovSubspaceDotProducts(preconditioner, inf.R, inf.R);
			inf.PP = CreateKrylovSubspaceDotProducts(preconditioner, inf.P, inf.P);
			inf.RP = CreateKrylovSubspaceDotProducts(preconditioner, inf.R, inf.P);
			return inf;
		}

		/** The coefficients of a vector's linear combination.
		 * They are the coefficients only from a linear combination.
		 * The result vector is a linear combination of vectors of 2 Krylov subspaces R and P.
		 * Coefficients of Krylov subspace R is in r list and for Krylov subspace P in p list.
		 * Lists can have fewer (even zero) elements than Krylov subspace vectors. */
		private class KSC
		{
			public KSC(int n)
			{
				r = new List<double>(n);
				p = new List<double>(n);
			}

			public List<double> r;
			public List<double> p;

			private static void Xupdate(List<double> x, double a, List<double> p)
			{
				int to2 = p.Count, to1 = Math.Min(to2, x.Count);
				for (int i = 0; i < to1; ++i)
					x[i] += a * p[i];
				for (int i = to1; i < to2; ++i)
					x.Add(a * p[i]);
			}

			/** It is the operation x += a * p.
			 * x and p are coefficients of Krylov subspace vectors linear combination. */
			public void Xupdate(double a, KSC p)
			{
				Xupdate(this.r, a, p.r);
				Xupdate(this.p, a, p.p);
			}

			private static void Rupdate(List<double> r, double a, List<double> p)
			{
				int to2 = p.Count, to1 = Math.Min(to2, r.Count - 1);
				for (int i = 0; i < to1; ++i)
					r[i + 1] -= a * p[i];
				if (to1 == -1) { r.Add(0); ++to1; }
				for (int i = to1; i < to2; ++i)
					r.Add(-a * p[i]);
			}

			/** It is the operation r -= a * A * p.
			 * x and p are coefficients of Krylov subspace vectors linear combination. */
			public void Rupdate(double a, KSC p)
			{
				Rupdate(this.r, a, p.r);
				Rupdate(this.p, a, p.p);
			}


			private static void Pupdate(List<double> r, double a, List<double> p)
			{
				for (int i = 0; i < p.Count; ++i) p[i] *= a;
				int to2 = r.Count, to1 = Math.Min(to2, p.Count);
				for (int i = 0; i < to1; ++i)
					p[i] += r[i];
				for (int i = to1; i < to2; ++i)
					p.Add(r[i]);
			}

			/** It is the operation p = r + a * p.
			 * r and p are coefficients of Krylov subspace vectors linear combination. */
			public void Pupdate(KSC r, double a)
			{
				Pupdate(r.r, a, this.r);
				Pupdate(r.p, a, this.p);
			}


			/** Calculates the linear combination of a vector.
			 * The operation is x = p_coefficints * P + r_coefficients * R
			 * <param name="inf">Block information with Krylov subspaces</param>
			 * <param name="xc">Coefficients of Krylov subspaces vectors</param>
			 * <return>The result vector</return>
			 */
			public IVector EvaluateVector(BlockInfo inf)
			{
				IVector x = inf.R[0].CreateZeroVectorWithSameFormat();
				for (int i = 0; i < r.Count; ++i)
					x.AxpyIntoThis(inf.R[i], r[i]);
				for (int i = 0; i < p.Count; ++i)
					x.AxpyIntoThis(inf.P[i], p[i]);
				return x;
			}

			private static double square1(List<double> RR, List<double> r)
			{
				double a = 0;
				for (int i = 0; i < r.Count; ++i)
					a += r[i] * r[i] * RR[2 * i];
				return a;
			}

			private static double square2(List<double> RR, List<double> r)
			{
				double a = 0;
				for (int i = 0; i < r.Count - 1; ++i)
					for (int j = i + 1; j < r.Count; j++)
						a += 2 * r[i] * r[j] * RR[i + j];
				return a;
			}

			private static double square3(List<double> RP, List<double> r, List<double> p)
			{
				double a = 0;
				for (int i = 0; i < r.Count; ++i)
					for (int j = 0; j < p.Count; j++)
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
			public double square(BlockInfo inf)
			{
				return square1(inf.RR, r) + square1(inf.PP, p) + square2(inf.RR, r) + square2(inf.PP, p) + square3(inf.RP, r, p);
			}


			private static double sandwich1(List<double> RR, List<double> r)
			{
				double a = 0;
				for (int i = 0; i < r.Count; ++i)
					a += r[i] * r[i] * RR[2 * i + 1];
				return a;
			}

			private static double sandwich2(List<double> RR, List<double> r)
			{
				double a = 0;
				for (int i = 0; i < r.Count - 1; ++i)
					for (int j = i + 1; j < r.Count; j++)
						a += 2 * r[i] * r[j] * RR[i + j + 1];
				return a;
			}

			private static double sandwich3(List<double> RP, List<double> r, List<double> p)
			{
				double a = 0;
				for (int i = 0; i < r.Count; ++i)
					for (int j = 0; j < p.Count; j++)
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
			public double sandwich(BlockInfo inf)
			{
				return sandwich1(inf.RR, r) + sandwich1(inf.PP, p) + sandwich2(inf.RR, r) + sandwich2(inf.PP, p) + sandwich3(inf.RP, r, p);
			}
		};



		private static readonly String AlgorithmName = "Block Proconditioned Conjugent Gradient";

		/** Extends iterative statistics. */
		public class IterativeStatistics : MGroup.LinearAlgebra.Iterative.IterativeStatistics
		{
			/// The square residual vector r^2
			public double resDotRes { get; set; }

			///  The solution of the algorithm
			public IVector solution { get; set; }
		}

		private static IterativeStatistics SolveSuccess(int iterations, double resDotRes, double resDotRes0, IVector solution)
		{
			return new IterativeStatistics
			{
				HasConverged = true,
				HasStagnated = false,
				AlgorithmName = AlgorithmName,
				NumIterationsRequired = iterations,
				ResidualNormRatioEstimation = resDotRes / resDotRes0,
				resDotRes = resDotRes,
				solution = solution,
			};
		}
		private static IterativeStatistics SolveFailure(int iterations, double resDotRes, double resDotRes0, IVector solution, bool stagnated)
		{
			return new IterativeStatistics
			{
				HasConverged = false,
				HasStagnated = stagnated,
				AlgorithmName = AlgorithmName,
				NumIterationsRequired = iterations,
				ResidualNormRatioEstimation = resDotRes / resDotRes0,
				resDotRes = resDotRes,
				solution = solution,
			};
		}


		/** Block Preconditioned Conjugent Gradient.
		 * It solves the problem A * x = b for x.
		 * <param name="A">The matrix A</param>
		 * <param name="preconditioner">The preconditioner. It can be null. In that case Jacobi preconditioner is used.</param>
		 * <param name="B">The vector b</param>
		 * <param name="X">The initial guess of solution. If null it is assumed zero.
		 * If not null, after algorithm's end, it contains the solution.</param>
		 * <param name="blockSize">Size of block. Usually 2 to 14 but above 8 starts to be unstable</param>
		 * <param name="maxIterationsProvider">Interface to stop algorithm because it exceeds maximum iterations</param>
		 * <param name="residualConvergence">Interface to stop algorithm because of successful convergence</param>
		 * <returns>Solution and statistics</returns>
		 */
		public static IterativeStatistics Solve(IMatrixView A, IPreconditioner preconditioner, IVectorView B, IVector X,
			int blockSize, IMaxIterationsProvider maxIterationsProvider, IResidualConvergence residualConvergence)
		{
			var linearTransformation = new ExplicitMatrixTransformation(A);
			if (preconditioner == null) preconditioner = new JacobiPreconditioner(A.GetDiagonalAsArray());
			return Solve(linearTransformation, preconditioner, B, X, blockSize, maxIterationsProvider, residualConvergence);
		}


		private static IVector EvaluateX(BlockInfo inf, IPreconditioner M, KSC x)
		{
			IVector X = inf.R[0].CreateZeroVectorWithSameFormat();
			M.SolveLinearSystem(x.EvaluateVector(inf), X);
			return X;
		}

		/** Block Preconditioned Conjugent Gradient.
		 * It solves the problem A * x = b for x.
		 * <param name="linearTransformation">A linear transformation matrix * vector corresponding to matrix A.</param>
		 * <param name="preconditioner">The preconditioner.</param>
		 * <param name="B">The vector b</param>
		 * <param name="X">The initial guess of solution. If null it is assumed zero.
		 * If not null, after algorithm's end, it contains the solution.</param>
		 * <param name="blockSize">Size of block. Usually 2 to 14 but above 8 starts to be unstable</param>
		 * <param name="maxIterationsProvider">Interface to stop algorithm because it exceeds maximum iterations</param>
		 * <param name="residualConvergence">Interface to stop algorithm because of successful convergence</param>
		 * <returns>Solution and statistics</returns>
		 */
		public static IterativeStatistics Solve(ILinearTransformation linearTransformation, IPreconditioner preconditioner, IVectorView B,
			IVector X, int blockSize, IMaxIterationsProvider maxIterationsProvider, IResidualConvergence residualConvergence)
		{
			Preconditions.CheckSystemSolutionDimensions(linearTransformation.NumRows, B.Length);
			Preconditions.CheckMultiplicationDimensions(linearTransformation.NumRows, linearTransformation.NumColumns);

			IVector R, P;
			if (X == null)
			{
				R = B.Copy();
				X = B.CreateZeroVectorWithSameFormat();
				Preconditions.CheckMultiplicationDimensions(linearTransformation.NumColumns, X.Length);
			}
			else
			{
				Preconditions.CheckMultiplicationDimensions(linearTransformation.NumColumns, X.Length);
				R = B.CreateZeroVectorWithSameFormat();
				linearTransformation.Multiply(X, R);
				R = B.Subtract(R);
			}
			BlockInfo inf = CreateBlockInfo(linearTransformation, preconditioner, R, blockSize);
			double rr = inf.RR[0], rr0 = rr;                          // rr = r * M * r
			if (residualConvergence.isConverged(rr)) return SolveSuccess(0, 1, 1, X);

			for (int step = 0; ;)
			{
				KSC r = new KSC(blockSize); r.r.Add(1);
				KSC p = new KSC(blockSize); p.p.Add(1);
				KSC x = new KSC(blockSize);

				for (int i = 0; i < blockSize; ++i)
				{
					double pAp = p.sandwich(inf);                     // pAp = p * M * A * M * p
					if (pAp <= 0)
					{
						if (i != 0) X.AddIntoThis(EvaluateX(inf, preconditioner, x));  // x = x0 + M * x. It multiplied with M, because it should be
						return SolveFailure(step * blockSize + i, rr, rr0, X, true);
					}
					double a = rr / pAp;
					x.Xupdate(a, p);                                  // x += a * p, but in should be x += a * M * p
					r.Rupdate(a, p);                                  // r -= a * A * p, but in should be r -= a * A * M * p
					double rr2 = r.square(inf);                       // rr2 = r * M * r
					double b = rr2 / rr;
					p.Pupdate(r, b);                                  // p = r + b * p, but in should be p = r + b * M * p
					rr = rr2;
					if (residualConvergence.isConverged(rr))
					{
						if (i != 0) X.AddIntoThis(EvaluateX(inf, preconditioner, x));  // x = x0 + M * x. It multiplied with M, because it should be
						return SolveSuccess(step * blockSize + i, rr, rr0, X);
					}
				}
				X.AddIntoThis(EvaluateX(inf, preconditioner, x));                 // x = x0 + M * x. It multiplied with M, because it should be
				++step;
				if (step * blockSize >= maxIterationsProvider.GetMaxIterations(linearTransformation.NumColumns))
					return SolveFailure(step * blockSize, rr, rr0, X, true);
				R = r.EvaluateVector(inf);                           // It didn't multiplied with M, because it shouldn't be
				P = p.EvaluateVector(inf);                           // It didn't multiplied with M, but it should be
				inf = CreateBlockInfo(linearTransformation, preconditioner, R, P, blockSize);
			}
		}




		private readonly IMaxIterationsProvider maxIterationsProvider = new PercentageMaxIterationsProvider(1);
		private readonly IResidualConvergence residualConvergence = new RegularResidualConvergence();
		private int blockSize;
		private double resDotRes;
		private IVector solution;

		/** Initialize some inner algorithm members.
		 * Use algorithm as an object only if static Solve() method does not fit on your needs.
		 * <param name="blockSize">Size of block. Usually 2 to 14 but above 8 starts to be unstable</param>
		 * <param name="maxIterationsProvider">Interface to stop algorithm because it exceeds maximum iterations</param>
		 * <param name="residualConvergence">Interface to stop algorithm because of successful convergence</param>
		 */
		public BPCCGAlgorithm(IMaxIterationsProvider maxIterationsProvider, IResidualConvergence residualConvergence, int blockSize = 6)
		{
			this.blockSize = blockSize;
			this.maxIterationsProvider = maxIterationsProvider;
			this.residualConvergence = residualConvergence;
		}

		/** Initialize some inner algorithm members.
		 * Use algorithm as an object only if static Solve() method does not fit on your needs.
		 * <param name="blockSize">Size of block. Usually 2 to 14 but above 8 starts to be unstable</param>
		 * <param name="residualConvergence">Interface to stop algorithm because of successful convergence</param>
		 */
		public BPCCGAlgorithm(IResidualConvergence residualConvergence, int blockSize = 6)
		{
			this.blockSize = blockSize;
			this.residualConvergence = residualConvergence;
		}

		/** Initialize some inner algorithm members.
		 * Use algorithm as an object only if static Solve() method does not fit on your needs.
		 * <param name="blockSize">Size of block. Usually 2 to 14 but above 8 starts to be unstable</param>
		 * <param name="maxIterationsProvider">Interface to stop algorithm because it exceeds maximum iterations</param>
		 */
		public BPCCGAlgorithm(IMaxIterationsProvider maxIterationsProvider, int blockSize = 6)
		{
			this.blockSize = blockSize;
			this.maxIterationsProvider = maxIterationsProvider;
		}

		/** Initialize some inner algorithm members.
		 * Use algorithm as an object only if static Solve() method does not fit on your needs.
		 * <param name="blockSize">Size of block. Usually 2 to 14 but above 8 starts to be unstable</param>
		 */
		public BPCCGAlgorithm(int blockSize = 6) { this.blockSize = blockSize; }


		/** Block Preconditioned Conjugent Gradient.
		 * It solves the problem A * x = b for x.
		 * <param name="A">The matrix A</param>
		 * <param name="preconditioner">The preconditioner. It can be null. In that case Jacobi preconditioner is used.</param>
		 * <param name="B">The vector b</param>
		 * <param name="X">The initial guess of solution. If null it is assumed zero.</param>
		 * <returns>Solution and statistics</returns>
		 */
		public IterativeStatistics Solve(IMatrixView A, IVectorView B, IPreconditioner preconditioner = null, IVectorView X = null)
		{
			IterativeStatistics stat = Solve(A, preconditioner, B, X == null ? null : X.Copy(), blockSize, maxIterationsProvider, residualConvergence);
			resDotRes = stat.resDotRes;
			solution = stat.solution;
			return stat;
		}

		/** Block Preconditioned Conjugent Gradient.
		 * It solves the problem A * x = b for x.
		 * <param name="linearTransformation">A linear transformation of A * M * vector where M is the inverse preconditioner matrix</param>
		 * <param name="preconditioner">The preconditioner.</param>
		 * <param name="B">The vector b</param>
		 * <param name="X">The initial guess of solution. If null it is assumed zero.</param>
		 * <returns>Solution and statistics</returns>
		 */
		public IterativeStatistics Solve(ILinearTransformation linearTransformation, IVectorView B, IPreconditioner preconditioner, IVectorView X = null)
		{
			IterativeStatistics stat = Solve(linearTransformation, preconditioner, B, X == null ? null : X.Copy(), blockSize, maxIterationsProvider, residualConvergence);
			resDotRes = stat.resDotRes;
			solution = stat.solution;
			return stat;
		}

		/** Block Preconditioned Conjugent Gradient.
		 * It solves the problem A * x = b for x.
		 * <param name="A">The matrix</param>
		 * <param name="M">The inverse preconditioner matrix</param>
		 * <param name="B">The vector b</param>
		 * <param name="X">The initial guess of solution. If null it is assumed zero.</param>
		 * <returns>Solution and statistics</returns>
		 */
		public IterativeStatistics Solve(IMatrixView A, IVectorView B, IVectorView X) => Solve(A, B, null, X);

		/// <summary>
		/// The dot product <see cref="Residual"/> * <see cref="Residual"/>.
		/// </summary>
		public double ResDotRes => resDotRes;

		/// <summary>
		/// The current approximation to the solution of the linear system A * x = b
		/// </summary>
		public IVectorView Solution => solution;

		public IMaxIterationsProvider MaxIterationsProvider => maxIterationsProvider;
		public IResidualConvergence ResidualConvergence => residualConvergence;
		public int BlockSize => blockSize;


		/// <summary>
		/// Releases references to the vectors and matrices used by this object and sets scalars to their default values.
		/// </summary>
		public void Clear() { solution = null; }
	}
}
