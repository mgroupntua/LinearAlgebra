#region Translated by Jose Antonio De Santiago-Castillo.

//Translated by Jose Antonio De Santiago-Castillo.
//E-mail:JAntonioDeSantiago@gmail.com
//Website: www.DotNumerics.com
//
//Fortran to C# Translation.
//Translated by:
//F2CSharp Version 0.72 (Dicember 7, 2009)
//Code Optimizations: , assignment operator, for-loop: array indexes
//
#endregion

using System;

namespace MGroup.LinearAlgebra.Implementations.DotNumerics
{
    /// <summary>
    /// -- LAPACK auxiliary routine (version 3.1) --
    /// Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
    /// November 2006
    /// Purpose
    /// =======
    /// 
    /// DLAQPS computes a step of QR factorization with column pivoting
    /// of a real M-by-N matrix A by using Blas-3.  It tries to factorize
    /// NB columns from A starting from the row OFFSET+1, and updates all
    /// of the matrix with Blas-3 xGEMM.
    /// 
    /// In some cases, due to catastrophic cancellations, it cannot
    /// factorize NB columns.  Hence, the actual number of factorized
    /// columns is returned in KB.
    /// 
    /// Block A(1:OFFSET,1:N) is accordingly pivoted, but not factorized.
    /// 
    ///</summary>
    public class DLAQPS
    {
    

        #region Dependencies
        
        DGEMM _dgemm; DGEMV _dgemv; DLARFG _dlarfg; DSWAP _dswap; IDAMAX _idamax; DLAMCH _dlamch; DNRM2 _dnrm2; 

        #endregion


        #region Variables
        
        const double ZERO = 0.0E+0; const double ONE = 1.0E+0; 

        #endregion

        public DLAQPS(DGEMM dgemm, DGEMV dgemv, DLARFG dlarfg, DSWAP dswap, IDAMAX idamax, DLAMCH dlamch, DNRM2 dnrm2)
        {
    

            #region Set Dependencies
            
            _dgemm = dgemm; _dgemv = dgemv; _dlarfg = dlarfg; _dswap = dswap; _idamax = idamax; 
            _dlamch = dlamch;_dnrm2 = dnrm2; 

            #endregion

        }
    
        public DLAQPS()
        {
    

            #region Dependencies (Initialization)
            
            var lsame = new LSAME();
            var xerbla = new XERBLA();
            var dlamc3 = new DLAMC3();
            var dlapy2 = new DLAPY2();
            var dnrm2 = new DNRM2();
            var dscal = new DSCAL();
            var dswap = new DSWAP();
            var idamax = new IDAMAX();
            var dgemm = new DGEMM(lsame, xerbla);
            var dgemv = new DGEMV(lsame, xerbla);
            var dlamc1 = new DLAMC1(dlamc3);
            var dlamc4 = new DLAMC4(dlamc3);
            var dlamc5 = new DLAMC5(dlamc3);
            var dlamc2 = new DLAMC2(dlamc3, dlamc1, dlamc4, dlamc5);
            var dlamch = new DLAMCH(lsame, dlamc2);
            var dlarfg = new DLARFG(dlamch, dlapy2, dnrm2, dscal);

            #endregion


            #region Set Dependencies
            
            _dgemm = dgemm; _dgemv = dgemv; _dlarfg = dlarfg; _dswap = dswap; _idamax = idamax; 
            _dlamch = dlamch;_dnrm2 = dnrm2; 

            #endregion

        }
        /// <summary>
        /// Purpose
        /// =======
        /// 
        /// DLAQPS computes a step of QR factorization with column pivoting
        /// of a real M-by-N matrix A by using Blas-3.  It tries to factorize
        /// NB columns from A starting from the row OFFSET+1, and updates all
        /// of the matrix with Blas-3 xGEMM.
        /// 
        /// In some cases, due to catastrophic cancellations, it cannot
        /// factorize NB columns.  Hence, the actual number of factorized
        /// columns is returned in KB.
        /// 
        /// Block A(1:OFFSET,1:N) is accordingly pivoted, but not factorized.
        /// 
        ///</summary>
        /// <param name="M">
        /// (input) INTEGER
        /// The number of rows of the matrix A. M .GE. 0.
        ///</param>
        /// <param name="N">
        /// (input) INTEGER
        /// The number of columns of the matrix A. N .GE. 0
        ///</param>
        /// <param name="OFFSET">
        /// (input) INTEGER
        /// The number of rows of A that have been factorized in
        /// previous steps.
        ///</param>
        /// <param name="NB">
        /// (input) INTEGER
        /// The number of columns to factorize.
        ///</param>
        /// <param name="KB">
        /// (output) INTEGER
        /// The number of columns actually factorized.
        ///</param>
        /// <param name="A">
        /// (input/output) DOUBLE PRECISION array, dimension (LDA,N)
        /// On entry, the M-by-N matrix A.
        /// On exit, block A(OFFSET+1:M,1:KB) is the triangular
        /// factor obtained and block A(1:OFFSET,1:N) has been
        /// accordingly pivoted, but no factorized.
        /// The rest of the matrix, block A(OFFSET+1:M,KB+1:N) has
        /// been updated.
        ///</param>
        /// <param name="LDA">
        /// (input) INTEGER
        /// The leading dimension of the array A. LDA .GE. max(1,M).
        ///</param>
        /// <param name="JPVT">
        /// (input/output) INTEGER array, dimension (N)
        /// JPVT(I) = K .LE.=.GT. Column K of the full matrix A has been
        /// permuted into position I in AP.
        ///</param>
        /// <param name="TAU">
        /// (output) DOUBLE PRECISION array, dimension (KB)
        /// The scalar factors of the elementary reflectors.
        ///</param>
        /// <param name="VN1">
        /// (input/output) DOUBLE PRECISION array, dimension (N)
        /// The vector with the partial column norms.
        ///</param>
        /// <param name="VN2">
        /// (input/output) DOUBLE PRECISION array, dimension (N)
        /// The vector with the exact column norms.
        ///</param>
        /// <param name="AUXV">
        /// (input/output) DOUBLE PRECISION array, dimension (NB)
        /// Auxiliar vector.
        ///</param>
        /// <param name="F">
        /// (input/output) DOUBLE PRECISION array, dimension (LDF,NB)
        /// Matrix F' = L*Y'*A.
        ///</param>
        /// <param name="LDF">
        /// (input) INTEGER
        /// The leading dimension of the array F. LDF .GE. max(1,N).
        ///</param>
        public void Run(int M, int N, int OFFSET, int NB, ref int KB, ref double[] A, int offset_a
                         , int LDA, ref int[] JPVT, int offset_jpvt, ref double[] TAU, int offset_tau, ref double[] VN1, int offset_vn1, ref double[] VN2, int offset_vn2, ref double[] AUXV, int offset_auxv
                         , ref double[] F, int offset_f, int LDF)
        {

            #region Variables
            
            var ITEMP = 0; var J = 0; var K = 0; var LASTRK = 0; var LSTICC = 0; var PVT = 0; var RK = 0; double AKK = 0; 
            double TEMP = 0;double TEMP2 = 0; double TOL3Z = 0; 

            #endregion


            #region Implicit Variables
            
            var F_K = 0; 

            #endregion


            #region Array Index Correction
            
             var o_a = -1 - LDA + offset_a;  var o_jpvt = -1 + offset_jpvt;  var o_tau = -1 + offset_tau; 
             var o_vn1 = -1 + offset_vn1; var o_vn2 = -1 + offset_vn2;  var o_auxv = -1 + offset_auxv; 
             var o_f = -1 - LDF + offset_f;

            #endregion


            #region Prolog
            
            // *
            // *  -- LAPACK auxiliary routine (version 3.1) --
            // *     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
            // *     November 2006
            // *
            // *     .. Scalar Arguments ..
            // *     ..
            // *     .. Array Arguments ..
            // *     ..
            // *
            // *  Purpose
            // *  =======
            // *
            // *  DLAQPS computes a step of QR factorization with column pivoting
            // *  of a real M-by-N matrix A by using Blas-3.  It tries to factorize
            // *  NB columns from A starting from the row OFFSET+1, and updates all
            // *  of the matrix with Blas-3 xGEMM.
            // *
            // *  In some cases, due to catastrophic cancellations, it cannot
            // *  factorize NB columns.  Hence, the actual number of factorized
            // *  columns is returned in KB.
            // *
            // *  Block A(1:OFFSET,1:N) is accordingly pivoted, but not factorized.
            // *
            // *  Arguments
            // *  =========
            // *
            // *  M       (input) INTEGER
            // *          The number of rows of the matrix A. M >= 0.
            // *
            // *  N       (input) INTEGER
            // *          The number of columns of the matrix A. N >= 0
            // *
            // *  OFFSET  (input) INTEGER
            // *          The number of rows of A that have been factorized in
            // *          previous steps.
            // *
            // *  NB      (input) INTEGER
            // *          The number of columns to factorize.
            // *
            // *  KB      (output) INTEGER
            // *          The number of columns actually factorized.
            // *
            // *  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
            // *          On entry, the M-by-N matrix A.
            // *          On exit, block A(OFFSET+1:M,1:KB) is the triangular
            // *          factor obtained and block A(1:OFFSET,1:N) has been
            // *          accordingly pivoted, but no factorized.
            // *          The rest of the matrix, block A(OFFSET+1:M,KB+1:N) has
            // *          been updated.
            // *
            // *  LDA     (input) INTEGER
            // *          The leading dimension of the array A. LDA >= max(1,M).
            // *
            // *  JPVT    (input/output) INTEGER array, dimension (N)
            // *          JPVT(I) = K <==> Column K of the full matrix A has been
            // *          permuted into position I in AP.
            // *
            // *  TAU     (output) DOUBLE PRECISION array, dimension (KB)
            // *          The scalar factors of the elementary reflectors.
            // *
            // *  VN1     (input/output) DOUBLE PRECISION array, dimension (N)
            // *          The vector with the partial column norms.
            // *
            // *  VN2     (input/output) DOUBLE PRECISION array, dimension (N)
            // *          The vector with the exact column norms.
            // *
            // *  AUXV    (input/output) DOUBLE PRECISION array, dimension (NB)
            // *          Auxiliar vector.
            // *
            // *  F       (input/output) DOUBLE PRECISION array, dimension (LDF,NB)
            // *          Matrix F' = L*Y'*A.
            // *
            // *  LDF     (input) INTEGER
            // *          The leading dimension of the array F. LDF >= max(1,N).
            // *
            // *  Further Details
            // *  ===============
            // *
            // *  Based on contributions by
            // *    G. Quintana-Orti, Depto. de Informatica, Universidad Jaime I, Spain
            // *    X. Sun, Computer Science Dept., Duke University, USA
            // *
            // *  Partial column norm updating strategy modified by
            // *    Z. Drmac and Z. Bujanovic, Dept. of Mathematics,
            // *    University of Zagreb, Croatia.
            // *    June 2006.
            // *  For more details see LAPACK Working Note 176.
            // *  =====================================================================
            // *
            // *     .. Parameters ..
            // *     ..
            // *     .. Local Scalars ..
            // *     ..
            // *     .. External Subroutines ..
            // *     ..
            // *     .. Intrinsic Functions ..
            //      INTRINSIC          ABS, DBLE, MAX, MIN, NINT, SQRT;
            // *     ..
            // *     .. External Functions ..
            // *     ..
            // *     .. Executable Statements ..
            // *

            #endregion


            #region Body
            
            LASTRK = Math.Min(M, N + OFFSET);
            LSTICC = 0;
            K = 0;
            TOL3Z = Math.Sqrt(_dlamch.Run("Epsilon"));
            // *
            // *     Beginning of while loop.
            // *
        LABEL10:;
            if (K < NB && LSTICC == 0)
            {
                K += 1;
                RK = OFFSET + K;
                // *
                // *        Determine ith pivot column and swap if necessary
                // *
                PVT = K - 1 + _idamax.Run(N - K + 1, VN1, K + o_vn1, 1);
                if (PVT != K)
                {
                    _dswap.Run(M, ref A, 1+PVT * LDA + o_a, 1, ref A, 1+K * LDA + o_a, 1);
                    _dswap.Run(K - 1, ref F, PVT+1 * LDF + o_f, LDF, ref F, K+1 * LDF + o_f, LDF);
                    ITEMP = JPVT[PVT + o_jpvt];
                    JPVT[PVT + o_jpvt] = JPVT[K + o_jpvt];
                    JPVT[K + o_jpvt] = ITEMP;
                    VN1[PVT + o_vn1] = VN1[K + o_vn1];
                    VN2[PVT + o_vn2] = VN2[K + o_vn2];
                }
                // *
                // *        Apply previous Householder reflectors to column K:
                // *        A(RK:M,K) := A(RK:M,K) - A(RK:M,1:K-1)*F(K,1:K-1)'.
                // *
                if (K > 1)
                {
                    _dgemv.Run("No transpose", M - RK + 1, K - 1,  - ONE, A, RK+1 * LDA + o_a, LDA
                                    , F, K+1 * LDF + o_f, LDF, ONE, ref A, RK+K * LDA + o_a, 1);
                }
                // *
                // *        Generate elementary reflector H(k).
                // *
                if (RK < M)
                {
                    _dlarfg.Run(M - RK + 1, ref A[RK+K * LDA + o_a], ref A, RK + 1+K * LDA + o_a, 1, ref TAU[K + o_tau]);
                }
                else
                {
                    _dlarfg.Run(1, ref A[RK+K * LDA + o_a], ref A, RK+K * LDA + o_a, 1, ref TAU[K + o_tau]);
                }
                // *
                AKK = A[RK+K * LDA + o_a];
                A[RK+K * LDA + o_a] = ONE;
                // *
                // *        Compute Kth column of F:
                // *
                // *        Compute  F(K+1:N,K) := tau(K)*A(RK:M,K+1:N)'*A(RK:M,K).
                // *
                if (K < N)
                {
                    _dgemv.Run("Transpose", M - RK + 1, N - K, TAU[K + o_tau], A, RK+(K + 1) * LDA + o_a, LDA
                                    , A, RK+K * LDA + o_a, 1, ZERO, ref F, K + 1+K * LDF + o_f, 1);
                }
                // *
                // *        Padding F(1:K,K) with zeros.
                // *
                F_K = K * LDF + o_f;
                for (J = 1; J <= K; J++)
                {
                    F[J + F_K] = ZERO;
                }
                // *
                // *        Incremental updating of F:
                // *        F(1:N,K) := F(1:N,K) - tau(K)*F(1:N,1:K-1)*A(RK:M,1:K-1)'
                // *                    *A(RK:M,K).
                // *
                if (K > 1)
                {
                    _dgemv.Run("Transpose", M - RK + 1, K - 1,  - TAU[K + o_tau], A, RK+1 * LDA + o_a, LDA
                                    , A, RK+K * LDA + o_a, 1, ZERO, ref AUXV, 1 + o_auxv, 1);
                    // *
                    _dgemv.Run("No transpose", N, K - 1, ONE, F, 1+1 * LDF + o_f, LDF
                                    , AUXV, 1 + o_auxv, 1, ONE, ref F, 1+K * LDF + o_f, 1);
                }
                // *
                // *        Update the current row of A:
                // *        A(RK,K+1:N) := A(RK,K+1:N) - A(RK,1:K)*F(K+1:N,1:K)'.
                // *
                if (K < N)
                {
                    _dgemv.Run("No transpose", N - K, K,  - ONE, F, K + 1+1 * LDF + o_f, LDF
                                    , A, RK+1 * LDA + o_a, LDA, ONE, ref A, RK+(K + 1) * LDA + o_a, LDA);
                }
                // *
                // *        Update partial column norms.
                // *
                if (RK < LASTRK)
                {
                    for (J = K + 1; J <= N; J++)
                    {
                        if (VN1[J + o_vn1] != ZERO)
                        {
                            // *
                            // *                 NOTE: The following 4 lines follow from the analysis in
                            // *                 Lapack Working Note 176.
                            // *
                            TEMP = Math.Abs(A[RK+J * LDA + o_a]) / VN1[J + o_vn1];
                            TEMP = Math.Max(ZERO, (ONE + TEMP) * (ONE - TEMP));
                            TEMP2 = TEMP * Math.Pow(VN1[J + o_vn1] / VN2[J + o_vn2],2);
                            if (TEMP2 <= TOL3Z)
                            {
                                VN2[J + o_vn2] = Convert.ToDouble(LSTICC);
                                LSTICC = J;
                            }
                            else
                            {
                                VN1[J + o_vn1] *= Math.Sqrt(TEMP);
                            }
                        }
                    }
                }
                // *
                A[RK+K * LDA + o_a] = AKK;
                // *
                // *        End of while loop.
                // *
                goto LABEL10;
            }
            KB = K;
            RK = OFFSET + KB;
            // *
            // *     Apply the block reflector to the rest of the matrix:
            // *     A(OFFSET+KB+1:M,KB+1:N) := A(OFFSET+KB+1:M,KB+1:N) -
            // *                         A(OFFSET+KB+1:M,1:KB)*F(KB+1:N,1:KB)'.
            // *
            if (KB < Math.Min(N, M - OFFSET))
            {
                _dgemm.Run("No transpose", "Transpose", M - RK, N - KB, KB,  - ONE
                                , A, RK + 1+1 * LDA + o_a, LDA, F, KB + 1+1 * LDF + o_f, LDF, ONE, ref A, RK + 1+(KB + 1) * LDA + o_a
                                , LDA);
            }
            // *
            // *     Recomputation of difficult columns.
            // *
        LABEL40:;
            if (LSTICC > 0)
            {
                ITEMP = (int)Math.Round(VN2[LSTICC + o_vn2]);
                VN1[LSTICC + o_vn1] = _dnrm2.Run(M - RK, A, RK + 1+LSTICC * LDA + o_a, 1);
                // *
                // *        NOTE: The computation of VN1( LSTICC ) relies on the fact that 
                // *        SNRM2 does not fail on vectors with norm below the value of
                // *        SQRT(DLAMCH('S')) 
                // *
                VN2[LSTICC + o_vn2] = VN1[LSTICC + o_vn1];
                LSTICC = ITEMP;
                goto LABEL40;
            }
            // *
            return;
            // *
            // *     End of DLAQPS
            // *

            #endregion

        }
    }
}
