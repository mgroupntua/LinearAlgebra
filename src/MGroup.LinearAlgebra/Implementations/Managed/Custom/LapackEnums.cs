namespace MGroup.LinearAlgebra.Implementations.Managed.Custom
{
	using System;

	public enum LapackDataType
	{
		SingleReal,    // S
		DoubleReal,    // D
		SingleComplex, // C
		DoubleComplex, // Z
	}

	public enum LapackMatrixType
	{
		General,      // GE
		Symmetric,    // SY
		Hermitian,    // HE
		PositiveDefinite, // PO
		BandGeneral,  // GB
		BandPositiveDefinite, // PB
		Triangular,   // TR
		Orthogonal,   // OR
		Unitary,      // UN
		Other,        // fallback
	}

	public enum LapackOperation
	{
		TRF,
		QRF,
		RQF,
		LQF,
		QLF,
		HRD,
		BRD,
		TRI,
		GST,
		UUM,
		EBZ,
		TRD,
		Other,
	}

	[Flags]
	public enum LapackOptions
	{
		None = 0,
		Upper = 1 << 0,
		Lower = 1 << 1,
		Transpose = 1 << 2,
		NoTranspose = 1 << 3,
		ConjugateTranspose = 1 << 4,
		UnitDiagonal = 1 << 5,
		NonUnitDiagonal = 1 << 6
	}

	public enum IlaenvSpec
	{
		OptimalBlockSize = 1,
		MinimumBlockSize = 2,
		CrossoverPoint = 3,
		NumberOfShifts = 4,
		MinimumColumnDimension = 5,
		SvdCrossover = 6,
		NumberOfProcessors = 7,
		MultishiftCrossover = 8,
		MaxSubproblemSize = 9,
		IeeeNaNCheck = 10,
		IeeeInfinityCheck = 11,
		IpArmq = 12, // 12–16 delegated
	}

	public enum DlamchParam
	{
		Epsilon,           // "E"
		SafeMinimum,       // "S"
		Base,              // "B"
		Precision,         // "P"
		MantissaDigits,    // "N"
		Rounding,          // "R"
		MinExponent,       // "M"
		UnderflowThreshold,// "U"
		MaxExponent,       // "L"
		OverflowThreshold  // "O"
	}
}
