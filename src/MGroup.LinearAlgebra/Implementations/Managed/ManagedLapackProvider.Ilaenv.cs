namespace MGroup.LinearAlgebra.Implementations.Managed
{
	using System;
	using System.Diagnostics;

	public partial class ManagedLapackProvider : ILapackProvider
	{
		/// <summary>
		/// Determines block size.
		/// </summary>
		private int Ilaenv(IlaenvSpec spec, LapackDataType type, LapackMatrixType matrix, LapackOperation op, LapackOptions opts, int n1, int n2, int n3, int n4)
		{
			bool isReal = type == LapackDataType.SingleReal || type == LapackDataType.DoubleReal;
			bool isComplex = !isReal;

			switch (spec)
			{
				// ============================================================
				// ISPEC = 1 → Optimal block size
				// ============================================================
				case IlaenvSpec.OptimalBlockSize:
					{
						switch (matrix)
						{
							case LapackMatrixType.General:
								switch (op)
								{
									case LapackOperation.TRF: return 64;
									case LapackOperation.QRF:
									case LapackOperation.RQF:
									case LapackOperation.LQF:
									case LapackOperation.QLF: return 32;
									case LapackOperation.HRD:
									case LapackOperation.BRD: return 32;
									case LapackOperation.TRI: return 64;
								}
								break;

							case LapackMatrixType.PositiveDefinite:
								if (op == LapackOperation.TRF) return 64;
								break;

							case LapackMatrixType.Symmetric:
								if (op == LapackOperation.TRF) return 64;
								if (isReal && op == LapackOperation.TRD) return 32;
								if (isReal && op == LapackOperation.GST) return 64;
								break;

							case LapackMatrixType.Hermitian:
								if (op == LapackOperation.TRF) return 64;
								if (op == LapackOperation.TRD) return 32;
								if (op == LapackOperation.GST) return 64;
								break;

							case LapackMatrixType.Orthogonal:
							case LapackMatrixType.Unitary:
								switch (op)
								{
									case LapackOperation.QRF:
									case LapackOperation.RQF:
									case LapackOperation.LQF:
									case LapackOperation.QLF:
									case LapackOperation.HRD:
									case LapackOperation.TRI:
									case LapackOperation.BRD:
										return 32;
								}
								break;

							case LapackMatrixType.BandGeneral:
								if (op == LapackOperation.TRF)
									return (n4 <= 64) ? 1 : 32;
								break;

							case LapackMatrixType.BandPositiveDefinite:
								if (op == LapackOperation.TRF)
									return (n2 <= 64) ? 1 : 32;
								break;

							case LapackMatrixType.Triangular:
								if (op == LapackOperation.TRI) return 64;
								break;

							case LapackMatrixType.Other:
								if (op == LapackOperation.UUM) return 64;
								if (op == LapackOperation.EBZ) return 1;
								break;
						}

						return 1;
					}

				// ============================================================
				// ISPEC = 2 → Minimum block size
				// ============================================================
				case IlaenvSpec.MinimumBlockSize:
					return 2;

				// ============================================================
				// ISPEC = 3 → Crossover point
				// ============================================================
				case IlaenvSpec.CrossoverPoint:
					{
						if (matrix == LapackMatrixType.General)
						{
							switch (op)
							{
								case LapackOperation.QRF:
								case LapackOperation.RQF:
								case LapackOperation.LQF:
								case LapackOperation.QLF:
								case LapackOperation.HRD:
								case LapackOperation.BRD:
									return 128;
							}
						}

						if ((matrix == LapackMatrixType.Symmetric || matrix == LapackMatrixType.Hermitian) && op == LapackOperation.TRD)
							return 32;

						if (matrix == LapackMatrixType.Orthogonal || matrix == LapackMatrixType.Unitary)
							return 128;

						return 0;
					}

				// ============================================================
				// ISPEC = 4–9 → constants
				// ============================================================
				case IlaenvSpec.NumberOfShifts: return 6;
				case IlaenvSpec.MinimumColumnDimension: return 2;
				case IlaenvSpec.SvdCrossover: return (int)(Math.Min(n1, n2) * 1.6);
				case IlaenvSpec.NumberOfProcessors: return 1;
				case IlaenvSpec.MultishiftCrossover: return 50;
				case IlaenvSpec.MaxSubproblemSize: return 25;

				// ============================================================
				// IEEE checks
				// ============================================================
				case IlaenvSpec.IeeeNaNCheck:
					return IeeeCheck(0);

				case IlaenvSpec.IeeeInfinityCheck:
					return IeeeCheck(1);

				// ============================================================
				// ISPEC = 12–16 → delegated to IPARMQ
				// ============================================================
				default:
					return Iparmq(spec, n1, n2, n3, n4);
			}
		}
	}
}
