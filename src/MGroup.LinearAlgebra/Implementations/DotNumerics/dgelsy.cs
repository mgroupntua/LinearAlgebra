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
    /// -- LAPACK driver routine (version 3.1) --
    /// Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
    /// November 2006
    /// Purpose
    /// =======
    /// 
    /// DGELSY computes the minimum-norm solution to a real linear least
    /// squares problem:
    /// minimize || A * X - B ||
    /// using a complete orthogonal factorization of A.  A is an M-by-N
    /// matrix which may be rank-deficient.
    /// 
    /// Several right hand side vectors b and solution vectors x can be
    /// handled in a single call; they are stored as the columns of the
    /// M-by-NRHS right hand side matrix B and the N-by-NRHS solution
    /// matrix X.
    /// 
    /// The routine first computes a QR factorization with column pivoting:
    /// A * P = Q * [ R11 R12 ]
    /// [  0  R22 ]
    /// with R11 defined as the largest leading submatrix whose estimated
    /// condition number is less than 1/RCOND.  The order of R11, RANK,
    /// is the effective rank of A.
    /// 
    /// Then, R22 is considered to be negligible, and R12 is annihilated
    /// by orthogonal transformations from the right, arriving at the
    /// complete orthogonal factorization:
    /// A * P = Q * [ T11 0 ] * Z
    /// [  0  0 ]
    /// The minimum-norm solution is then
    /// X = P * Z' [ inv(T11)*Q1'*B ]
    /// [        0       ]
    /// where Q1 consists of the first RANK columns of Q.
    /// 
    /// This routine is basically identical to the original xGELSX except
    /// three differences:
    /// o The call to the subroutine xGEQPF has been substituted by the
    /// the call to the subroutine xGEQP3. This subroutine is a Blas-3
    /// version of the QR factorization with column pivoting.
    /// o Matrix B (the right hand side) is updated with Blas-3.
    /// o The permutation of matrix B (the right hand side) is faster and
    /// more simple.
    /// 
    ///</summary>
    public class DGELSY
    {
    

        #region Dependencies
        
        ILAENV _ilaenv; DLAMCH _dlamch; DLANGE _dlange; DCOPY _dcopy; DGEQP3 _dgeqp3; DLABAD _dlabad; DLAIC1 _dlaic1; 
        DLASCL _dlascl;DLASET _dlaset; DORMQR _dormqr; DORMRZ _dormrz; DTRSM _dtrsm; DTZRZF _dtzrzf; XERBLA _xerbla; 

        #endregion


        #region Variables
        
        const int IMAX = 1; const int IMIN = 2; const double ZERO = 0.0E+0; const double ONE = 1.0E+0; 

        #endregion

        public DGELSY(ILAENV ilaenv, DLAMCH dlamch, DLANGE dlange, DCOPY dcopy, DGEQP3 dgeqp3, DLABAD dlabad, DLAIC1 dlaic1, DLASCL dlascl, DLASET dlaset, DORMQR dormqr
                      , DORMRZ dormrz, DTRSM dtrsm, DTZRZF dtzrzf, XERBLA xerbla)
        {
    

            #region Set Dependencies
            
            _ilaenv = ilaenv; _dlamch = dlamch; _dlange = dlange; _dcopy = dcopy; _dgeqp3 = dgeqp3; 
            _dlabad = dlabad;_dlaic1 = dlaic1; _dlascl = dlascl; _dlaset = dlaset; _dormqr = dormqr; 
            _dormrz = dormrz;_dtrsm = dtrsm; _dtzrzf = dtzrzf; _xerbla = xerbla; 

            #endregion

        }
    
        public DGELSY()
        {
    

            #region Dependencies (Initialization)
            
            var ieeeck = new IEEECK();
            var iparmq = new IPARMQ();
            var lsame = new LSAME();
            var dlamc3 = new DLAMC3();
            var dlassq = new DLASSQ();
            var dcopy = new DCOPY();
            var xerbla = new XERBLA();
            var dlapy2 = new DLAPY2();
            var dnrm2 = new DNRM2();
            var dscal = new DSCAL();
            var dswap = new DSWAP();
            var idamax = new IDAMAX();
            var dlabad = new DLABAD();
            var ddot = new DDOT();
            var daxpy = new DAXPY();
            var ilaenv = new ILAENV(ieeeck, iparmq);
            var dlamc1 = new DLAMC1(dlamc3);
            var dlamc4 = new DLAMC4(dlamc3);
            var dlamc5 = new DLAMC5(dlamc3);
            var dlamc2 = new DLAMC2(dlamc3, dlamc1, dlamc4, dlamc5);
            var dlamch = new DLAMCH(lsame, dlamc2);
            var dlange = new DLANGE(dlassq, lsame);
            var dgemv = new DGEMV(lsame, xerbla);
            var dger = new DGER(xerbla);
            var dlarf = new DLARF(dgemv, dger, lsame);
            var dlarfg = new DLARFG(dlamch, dlapy2, dnrm2, dscal);
            var dgeqr2 = new DGEQR2(dlarf, dlarfg, xerbla);
            var dgemm = new DGEMM(lsame, xerbla);
            var dtrmm = new DTRMM(lsame, xerbla);
            var dlarfb = new DLARFB(lsame, dcopy, dgemm, dtrmm);
            var dtrmv = new DTRMV(lsame, xerbla);
            var dlarft = new DLARFT(dgemv, dtrmv, lsame);
            var dgeqrf = new DGEQRF(dgeqr2, dlarfb, dlarft, xerbla, ilaenv);
            var dlaqp2 = new DLAQP2(dlarf, dlarfg, dswap, idamax, dlamch, dnrm2);
            var dlaqps = new DLAQPS(dgemm, dgemv, dlarfg, dswap, idamax, dlamch, dnrm2);
            var dorm2r = new DORM2R(lsame, dlarf, xerbla);
            var dormqr = new DORMQR(lsame, ilaenv, dlarfb, dlarft, dorm2r, xerbla);
            var dgeqp3 = new DGEQP3(dgeqrf, dlaqp2, dlaqps, dormqr, dswap, xerbla, ilaenv, dnrm2);
            var dlaic1 = new DLAIC1(ddot, dlamch);
            var dlascl = new DLASCL(lsame, dlamch, xerbla);
            var dlaset = new DLASET(lsame);
            var dlarzb = new DLARZB(lsame, dcopy, dgemm, dtrmm, xerbla);
            var dlarzt = new DLARZT(dgemv, dtrmv, xerbla, lsame);
            var dlarz = new DLARZ(daxpy, dcopy, dgemv, dger, lsame);
            var dormr3 = new DORMR3(lsame, dlarz, xerbla);
            var dormrz = new DORMRZ(lsame, ilaenv, dlarzb, dlarzt, dormr3, xerbla);
            var dtrsm = new DTRSM(lsame, xerbla);
            var dlatrz = new DLATRZ(dlarfg, dlarz);
            var dtzrzf = new DTZRZF(dlarzb, dlarzt, dlatrz, xerbla, ilaenv);

            #endregion


            #region Set Dependencies
            
            _ilaenv = ilaenv; _dlamch = dlamch; _dlange = dlange; _dcopy = dcopy; _dgeqp3 = dgeqp3; 
            _dlabad = dlabad;_dlaic1 = dlaic1; _dlascl = dlascl; _dlaset = dlaset; _dormqr = dormqr; 
            _dormrz = dormrz;_dtrsm = dtrsm; _dtzrzf = dtzrzf; _xerbla = xerbla; 

            #endregion

        }
        /// <summary>
        /// Purpose
        /// =======
        /// 
        /// DGELSY computes the minimum-norm solution to a real linear least
        /// squares problem:
        /// minimize || A * X - B ||
        /// using a complete orthogonal factorization of A.  A is an M-by-N
        /// matrix which may be rank-deficient.
        /// 
        /// Several right hand side vectors b and solution vectors x can be
        /// handled in a single call; they are stored as the columns of the
        /// M-by-NRHS right hand side matrix B and the N-by-NRHS solution
        /// matrix X.
        /// 
        /// The routine first computes a QR factorization with column pivoting:
        /// A * P = Q * [ R11 R12 ]
        /// [  0  R22 ]
        /// with R11 defined as the largest leading submatrix whose estimated
        /// condition number is less than 1/RCOND.  The order of R11, RANK,
        /// is the effective rank of A.
        /// 
        /// Then, R22 is considered to be negligible, and R12 is annihilated
        /// by orthogonal transformations from the right, arriving at the
        /// complete orthogonal factorization:
        /// A * P = Q * [ T11 0 ] * Z
        /// [  0  0 ]
        /// The minimum-norm solution is then
        /// X = P * Z' [ inv(T11)*Q1'*B ]
        /// [        0       ]
        /// where Q1 consists of the first RANK columns of Q.
        /// 
        /// This routine is basically identical to the original xGELSX except
        /// three differences:
        /// o The call to the subroutine xGEQPF has been substituted by the
        /// the call to the subroutine xGEQP3. This subroutine is a Blas-3
        /// version of the QR factorization with column pivoting.
        /// o Matrix B (the right hand side) is updated with Blas-3.
        /// o The permutation of matrix B (the right hand side) is faster and
        /// more simple.
        /// 
        ///</summary>
        /// <param name="M">
        /// (input) INTEGER
        /// The number of rows of the matrix A.  M .GE. 0.
        ///</param>
        /// <param name="N">
        /// (input) INTEGER
        /// The number of columns of the matrix A.  N .GE. 0.
        ///</param>
        /// <param name="NRHS">
        /// (input) INTEGER
        /// The number of right hand sides, i.e., the number of
        /// columns of matrices B and X. NRHS .GE. 0.
        ///</param>
        /// <param name="A">
        /// (input/output) DOUBLE PRECISION array, dimension (LDA,N)
        /// On entry, the M-by-N matrix A.
        /// On exit, A has been overwritten by details of its
        /// complete orthogonal factorization.
        ///</param>
        /// <param name="LDA">
        /// (input) INTEGER
        /// The leading dimension of the array A.  LDA .GE. max(1,M).
        ///</param>
        /// <param name="B">
        /// (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
        /// On entry, the M-by-NRHS right hand side matrix B.
        /// On exit, the N-by-NRHS solution matrix X.
        ///</param>
        /// <param name="LDB">
        /// (input) INTEGER
        /// The leading dimension of the array B. LDB .GE. max(1,M,N).
        ///</param>
        /// <param name="JPVT">
        /// (input/output) INTEGER array, dimension (N)
        /// On entry, if JPVT(i) .ne. 0, the i-th column of A is permuted
        /// to the front of AP, otherwise column i is a free column.
        /// On exit, if JPVT(i) = k, then the i-th column of AP
        /// was the k-th column of A.
        ///</param>
        /// <param name="RCOND">
        /// (input) DOUBLE PRECISION
        /// RCOND is used to determine the effective rank of A, which
        /// is defined as the order of the largest leading triangular
        /// submatrix R11 in the QR factorization with pivoting of A,
        /// whose estimated condition number .LT. 1/RCOND.
        ///</param>
        /// <param name="RANK">
        /// (output) INTEGER
        /// The effective rank of A, i.e., the order of the submatrix
        /// R11.  This is the same as the order of the submatrix T11
        /// in the complete orthogonal factorization of A.
        ///</param>
        /// <param name="WORK">
        /// (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
        /// On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
        ///</param>
        /// <param name="LWORK">
        /// (input) INTEGER
        /// The dimension of the array WORK.
        /// The unblocked strategy requires that:
        /// LWORK .GE. MAX( MN+3*N+1, 2*MN+NRHS ),
        /// where MN = min( M, N ).
        /// The block algorithm requires that:
        /// LWORK .GE. MAX( MN+2*N+NB*(N+1), 2*MN+NB*NRHS ),
        /// where NB is an upper bound on the blocksize returned
        /// by ILAENV for the routines DGEQP3, DTZRZF, STZRQF, DORMQR,
        /// and DORMRZ.
        /// 
        /// If LWORK = -1, then a workspace query is assumed; the routine
        /// only calculates the optimal size of the WORK array, returns
        /// this value as the first entry of the WORK array, and no error
        /// message related to LWORK is issued by XERBLA.
        ///</param>
        /// <param name="INFO">
        /// (output) INTEGER
        /// = 0: successful exit
        /// .LT. 0: If INFO = -i, the i-th argument had an illegal value.
        ///</param>
        public void Run(int M, int N, int NRHS, ref double[] A, int offset_a, int LDA, ref double[] B, int offset_b
                         , int LDB, ref int[] JPVT, int offset_jpvt, double RCOND, ref int RANK, ref double[] WORK, int offset_work, int LWORK
                         , ref int INFO)
        {

            #region Variables
            
            var LQUERY = false; var I = 0; var IASCL = 0; var IBSCL = 0; var ISMAX = 0; var ISMIN = 0; var J = 0; var LWKMIN = 0; 
            var LWKOPT = 0;var MN = 0; var NB = 0; var NB1 = 0; var NB2 = 0; var NB3 = 0; var NB4 = 0; double ANRM = 0; 
            double BIGNUM = 0;double BNRM = 0; double C1 = 0; double C2 = 0; double S1 = 0; double S2 = 0; double SMAX = 0; 
            double SMAXPR = 0;double SMIN = 0; double SMINPR = 0; double SMLNUM = 0; double WSIZE = 0; 

            #endregion


            #region Implicit Variables
            
            var B_J = 0; 

            #endregion


            #region Array Index Correction
            
             var o_a = -1 - LDA + offset_a;  var o_b = -1 - LDB + offset_b;  var o_jpvt = -1 + offset_jpvt; 
             var o_work = -1 + offset_work;

            #endregion


            #region Prolog
            
            // *
            // *  -- LAPACK driver routine (version 3.1) --
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
            // *  DGELSY computes the minimum-norm solution to a real linear least
            // *  squares problem:
            // *      minimize || A * X - B ||
            // *  using a complete orthogonal factorization of A.  A is an M-by-N
            // *  matrix which may be rank-deficient.
            // *
            // *  Several right hand side vectors b and solution vectors x can be
            // *  handled in a single call; they are stored as the columns of the
            // *  M-by-NRHS right hand side matrix B and the N-by-NRHS solution
            // *  matrix X.
            // *
            // *  The routine first computes a QR factorization with column pivoting:
            // *      A * P = Q * [ R11 R12 ]
            // *                  [  0  R22 ]
            // *  with R11 defined as the largest leading submatrix whose estimated
            // *  condition number is less than 1/RCOND.  The order of R11, RANK,
            // *  is the effective rank of A.
            // *
            // *  Then, R22 is considered to be negligible, and R12 is annihilated
            // *  by orthogonal transformations from the right, arriving at the
            // *  complete orthogonal factorization:
            // *     A * P = Q * [ T11 0 ] * Z
            // *                 [  0  0 ]
            // *  The minimum-norm solution is then
            // *     X = P * Z' [ inv(T11)*Q1'*B ]
            // *                [        0       ]
            // *  where Q1 consists of the first RANK columns of Q.
            // *
            // *  This routine is basically identical to the original xGELSX except
            // *  three differences:
            // *    o The call to the subroutine xGEQPF has been substituted by the
            // *      the call to the subroutine xGEQP3. This subroutine is a Blas-3
            // *      version of the QR factorization with column pivoting.
            // *    o Matrix B (the right hand side) is updated with Blas-3.
            // *    o The permutation of matrix B (the right hand side) is faster and
            // *      more simple.
            // *
            // *  Arguments
            // *  =========
            // *
            // *  M       (input) INTEGER
            // *          The number of rows of the matrix A.  M >= 0.
            // *
            // *  N       (input) INTEGER
            // *          The number of columns of the matrix A.  N >= 0.
            // *
            // *  NRHS    (input) INTEGER
            // *          The number of right hand sides, i.e., the number of
            // *          columns of matrices B and X. NRHS >= 0.
            // *
            // *  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
            // *          On entry, the M-by-N matrix A.
            // *          On exit, A has been overwritten by details of its
            // *          complete orthogonal factorization.
            // *
            // *  LDA     (input) INTEGER
            // *          The leading dimension of the array A.  LDA >= max(1,M).
            // *
            // *  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
            // *          On entry, the M-by-NRHS right hand side matrix B.
            // *          On exit, the N-by-NRHS solution matrix X.
            // *
            // *  LDB     (input) INTEGER
            // *          The leading dimension of the array B. LDB >= max(1,M,N).
            // *
            // *  JPVT    (input/output) INTEGER array, dimension (N)
            // *          On entry, if JPVT(i) .ne. 0, the i-th column of A is permuted
            // *          to the front of AP, otherwise column i is a free column.
            // *          On exit, if JPVT(i) = k, then the i-th column of AP
            // *          was the k-th column of A.
            // *
            // *  RCOND   (input) DOUBLE PRECISION
            // *          RCOND is used to determine the effective rank of A, which
            // *          is defined as the order of the largest leading triangular
            // *          submatrix R11 in the QR factorization with pivoting of A,
            // *          whose estimated condition number < 1/RCOND.
            // *
            // *  RANK    (output) INTEGER
            // *          The effective rank of A, i.e., the order of the submatrix
            // *          R11.  This is the same as the order of the submatrix T11
            // *          in the complete orthogonal factorization of A.
            // *
            // *  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
            // *          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
            // *
            // *  LWORK   (input) INTEGER
            // *          The dimension of the array WORK.
            // *          The unblocked strategy requires that:
            // *             LWORK >= MAX( MN+3*N+1, 2*MN+NRHS ),
            // *          where MN = min( M, N ).
            // *          The block algorithm requires that:
            // *             LWORK >= MAX( MN+2*N+NB*(N+1), 2*MN+NB*NRHS ),
            // *          where NB is an upper bound on the blocksize returned
            // *          by ILAENV for the routines DGEQP3, DTZRZF, STZRQF, DORMQR,
            // *          and DORMRZ.
            // *
            // *          If LWORK = -1, then a workspace query is assumed; the routine
            // *          only calculates the optimal size of the WORK array, returns
            // *          this value as the first entry of the WORK array, and no error
            // *          message related to LWORK is issued by XERBLA.
            // *
            // *  INFO    (output) INTEGER
            // *          = 0: successful exit
            // *          < 0: If INFO = -i, the i-th argument had an illegal value.
            // *
            // *  Further Details
            // *  ===============
            // *
            // *  Based on contributions by
            // *    A. Petitet, Computer Science Dept., Univ. of Tenn., Knoxville, USA
            // *    E. Quintana-Orti, Depto. de Informatica, Universidad Jaime I, Spain
            // *    G. Quintana-Orti, Depto. de Informatica, Universidad Jaime I, Spain
            // *
            // *  =====================================================================
            // *
            // *     .. Parameters ..
            // *     ..
            // *     .. Local Scalars ..
            // *     ..
            // *     .. External Functions ..
            // *     ..
            // *     .. External Subroutines ..
            // *     ..
            // *     .. Intrinsic Functions ..
            //      INTRINSIC          ABS, MAX, MIN;
            // *     ..
            // *     .. Executable Statements ..
            // *

            #endregion


            #region Body
            
            MN = Math.Min(M, N);
            ISMIN = MN + 1;
            ISMAX = 2 * MN + 1;
            // *
            // *     Test the input arguments.
            // *
            INFO = 0;
            LQUERY = LWORK ==  - 1;
            if (M < 0)
            {
                INFO =  - 1;
            }
            else
            {
                if (N < 0)
                {
                    INFO =  - 2;
                }
                else
                {
                    if (NRHS < 0)
                    {
                        INFO =  - 3;
                    }
                    else
                    {
                        if (LDA < Math.Max(1, M))
                        {
                            INFO =  - 5;
                        }
                        else
                        {
                            if (LDB < Math.Max(1, Math.Max(M, N)))
                            {
                                INFO =  - 7;
                            }
                        }
                    }
                }
            }
            // *
            // *     Figure out optimal block size
            // *
            if (INFO == 0)
            {
                if (MN == 0 || NRHS == 0)
                {
                    LWKMIN = 1;
                    LWKOPT = 1;
                }
                else
                {
                    NB1 = _ilaenv.Run(1, "DGEQRF", " ", M, N,  - 1,  - 1);
                    NB2 = _ilaenv.Run(1, "DGERQF", " ", M, N,  - 1,  - 1);
                    NB3 = _ilaenv.Run(1, "DORMQR", " ", M, N, NRHS,  - 1);
                    NB4 = _ilaenv.Run(1, "DORMRQ", " ", M, N, NRHS,  - 1);
                    NB = Math.Max(NB1, Math.Max(NB2, Math.Max(NB3, NB4)));
                    LWKMIN = MN + Math.Max(2 * MN, Math.Max(N + 1, MN + NRHS));
                    LWKOPT = Math.Max(LWKMIN, Math.Max(MN + 2 * N + NB * (N + 1), 2 * MN + NB * NRHS));
                }
                WORK[1 + o_work] = LWKOPT;
                // *
                if (LWORK < LWKMIN && !LQUERY)
                {
                    INFO =  - 12;
                }
            }
            // *
            if (INFO != 0)
            {
                _xerbla.Run("DGELSY",  - INFO);
                return;
            }
            else
            {
                if (LQUERY)
                {
                    return;
                }
            }
            // *
            // *     Quick return if possible
            // *
            if (MN == 0 || NRHS == 0)
            {
                RANK = 0;
                return;
            }
            // *
            // *     Get machine parameters
            // *
            SMLNUM = _dlamch.Run("S") / _dlamch.Run("P");
            BIGNUM = ONE / SMLNUM;
            _dlabad.Run(ref SMLNUM, ref BIGNUM);
            // *
            // *     Scale A, B if max entries outside range [SMLNUM,BIGNUM]
            // *
            ANRM = _dlange.Run("M", M, N, A, offset_a, LDA, ref WORK, offset_work);
            IASCL = 0;
            if (ANRM > ZERO && ANRM < SMLNUM)
            {
                // *
                // *        Scale matrix norm up to SMLNUM
                // *
                _dlascl.Run("G", 0, 0, ANRM, SMLNUM, M
                                 , N, ref A, offset_a, LDA, ref INFO);
                IASCL = 1;
            }
            else
            {
                if (ANRM > BIGNUM)
                {
                    // *
                    // *        Scale matrix norm down to BIGNUM
                    // *
                    _dlascl.Run("G", 0, 0, ANRM, BIGNUM, M
                                     , N, ref A, offset_a, LDA, ref INFO);
                    IASCL = 2;
                }
                else
                {
                    if (ANRM == ZERO)
                    {
                        // *
                        // *        Matrix all zero. Return zero solution.
                        // *
                        _dlaset.Run("F", Math.Max(M, N), NRHS, ZERO, ZERO, ref B, offset_b
                                         , LDB);
                        RANK = 0;
                        goto LABEL70;
                    }
                }
            }
            // *
            BNRM = _dlange.Run("M", M, NRHS, B, offset_b, LDB, ref WORK, offset_work);
            IBSCL = 0;
            if (BNRM > ZERO && BNRM < SMLNUM)
            {
                // *
                // *        Scale matrix norm up to SMLNUM
                // *
                _dlascl.Run("G", 0, 0, BNRM, SMLNUM, M
                                 , NRHS, ref B, offset_b, LDB, ref INFO);
                IBSCL = 1;
            }
            else
            {
                if (BNRM > BIGNUM)
                {
                    // *
                    // *        Scale matrix norm down to BIGNUM
                    // *
                    _dlascl.Run("G", 0, 0, BNRM, BIGNUM, M
                                     , NRHS, ref B, offset_b, LDB, ref INFO);
                    IBSCL = 2;
                }
            }
            // *
            // *     Compute QR factorization with column pivoting of A:
            // *        A * P = Q * R
            // *
            _dgeqp3.Run(M, N, ref A, offset_a, LDA, ref JPVT, offset_jpvt, ref WORK, 1 + o_work
                             , ref WORK, MN + 1 + o_work, LWORK - MN, ref INFO);
            WSIZE = MN + WORK[MN + 1 + o_work];
            // *
            // *     workspace: MN+2*N+NB*(N+1).
            // *     Details of Householder rotations stored in WORK(1:MN).
            // *
            // *     Determine RANK using incremental condition estimation
            // *
            WORK[ISMIN + o_work] = ONE;
            WORK[ISMAX + o_work] = ONE;
            SMAX = Math.Abs(A[1+1 * LDA + o_a]);
            SMIN = SMAX;
            if (Math.Abs(A[1+1 * LDA + o_a]) == ZERO)
            {
                RANK = 0;
                _dlaset.Run("F", Math.Max(M, N), NRHS, ZERO, ZERO, ref B, offset_b
                                 , LDB);
                goto LABEL70;
            }
            else
            {
                RANK = 1;
            }
            // *
        LABEL10:;
            if (RANK < MN)
            {
                I = RANK + 1;
                _dlaic1.Run(IMIN, RANK, WORK, ISMIN + o_work, SMIN, A, 1+I * LDA + o_a, A[I+I * LDA + o_a]
                                 , ref SMINPR, ref S1, ref C1);
                _dlaic1.Run(IMAX, RANK, WORK, ISMAX + o_work, SMAX, A, 1+I * LDA + o_a, A[I+I * LDA + o_a]
                                 , ref SMAXPR, ref S2, ref C2);
                // *
                if (SMAXPR * RCOND <= SMINPR)
                {
                    for (I = 1; I <= RANK; I++)
                    {
                        WORK[ISMIN + I - 1 + o_work] *= S1;
                        WORK[ISMAX + I - 1 + o_work] *= S2;
                    }
                    WORK[ISMIN + RANK + o_work] = C1;
                    WORK[ISMAX + RANK + o_work] = C2;
                    SMIN = SMINPR;
                    SMAX = SMAXPR;
                    RANK += 1;
                    goto LABEL10;
                }
            }
            // *
            // *     workspace: 3*MN.
            // *
            // *     Logically partition R = [ R11 R12 ]
            // *                             [  0  R22 ]
            // *     where R11 = R(1:RANK,1:RANK)
            // *
            // *     [R11,R12] = [ T11, 0 ] * Y
            // *
            if (RANK < N)
            {
                _dtzrzf.Run(RANK, N, ref A, offset_a, LDA, ref WORK, MN + 1 + o_work, ref WORK, 2 * MN + 1 + o_work
                                 , LWORK - 2 * MN, ref INFO);
            }
            // *
            // *     workspace: 2*MN.
            // *     Details of Householder rotations stored in WORK(MN+1:2*MN)
            // *
            // *     B(1:M,1:NRHS) := Q' * B(1:M,1:NRHS)
            // *
            _dormqr.Run("Left", "Transpose", M, NRHS, MN, ref A, offset_a
                             , LDA, WORK, 1 + o_work, ref B, offset_b, LDB, ref WORK, 2 * MN + 1 + o_work, LWORK - 2 * MN
                             , ref INFO);
            WSIZE = Math.Max(WSIZE, 2 * MN + WORK[2 * MN + 1 + o_work]);
            // *
            // *     workspace: 2*MN+NB*NRHS.
            // *
            // *     B(1:RANK,1:NRHS) := inv(T11) * B(1:RANK,1:NRHS)
            // *
            _dtrsm.Run("Left", "Upper", "No transpose", "Non-unit", RANK, NRHS
                            , ONE, A, offset_a, LDA, ref B, offset_b, LDB);
            // *
            for (J = 1; J <= NRHS; J++)
            {
                B_J = J * LDB + o_b;
                for (I = RANK + 1; I <= N; I++)
                {
                    B[I + B_J] = ZERO;
                }
            }
            // *
            // *     B(1:N,1:NRHS) := Y' * B(1:N,1:NRHS)
            // *
            if (RANK < N)
            {
                _dormrz.Run("Left", "Transpose", N, NRHS, RANK, N - RANK
                                 , A, offset_a, LDA, WORK, MN + 1 + o_work, ref B, offset_b, LDB, ref WORK, 2 * MN + 1 + o_work
                                 , LWORK - 2 * MN, ref INFO);
            }
            // *
            // *     workspace: 2*MN+NRHS.
            // *
            // *     B(1:N,1:NRHS) := P * B(1:N,1:NRHS)
            // *
            for (J = 1; J <= NRHS; J++)
            {
                B_J = J * LDB + o_b;
                for (I = 1; I <= N; I++)
                {
                    WORK[JPVT[I + o_jpvt] + o_work] = B[I + B_J];
                }
                _dcopy.Run(N, WORK, 1 + o_work, 1, ref B, 1+J * LDB + o_b, 1);
            }
            // *
            // *     workspace: N.
            // *
            // *     Undo scaling
            // *
            if (IASCL == 1)
            {
                _dlascl.Run("G", 0, 0, ANRM, SMLNUM, N
                                 , NRHS, ref B, offset_b, LDB, ref INFO);
                _dlascl.Run("U", 0, 0, SMLNUM, ANRM, RANK
                                 , RANK, ref A, offset_a, LDA, ref INFO);
            }
            else
            {
                if (IASCL == 2)
                {
                    _dlascl.Run("G", 0, 0, ANRM, BIGNUM, N
                                     , NRHS, ref B, offset_b, LDB, ref INFO);
                    _dlascl.Run("U", 0, 0, BIGNUM, ANRM, RANK
                                     , RANK, ref A, offset_a, LDA, ref INFO);
                }
            }
            if (IBSCL == 1)
            {
                _dlascl.Run("G", 0, 0, SMLNUM, BNRM, N
                                 , NRHS, ref B, offset_b, LDB, ref INFO);
            }
            else
            {
                if (IBSCL == 2)
                {
                    _dlascl.Run("G", 0, 0, BIGNUM, BNRM, N
                                     , NRHS, ref B, offset_b, LDB, ref INFO);
                }
            }
            // *
        LABEL70:;
            WORK[1 + o_work] = LWKOPT;
            // *
            return;
            // *
            // *     End of DGELSY
            // *

            #endregion

        }
    }
}
