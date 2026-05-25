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
    /// -- LAPACK routine (version 3.1) --
    /// Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
    /// November 2006
    /// Purpose
    /// =======
    /// 
    /// DGERQF computes an RQ factorization of a real M-by-N matrix A:
    /// A = R * Q.
    /// 
    ///</summary>
    public class DGERQF
    {
    

        #region Dependencies
        
        DGERQ2 _dgerq2; DLARFB _dlarfb; DLARFT _dlarft; XERBLA _xerbla; ILAENV _ilaenv; 

        #endregion

        public DGERQF(DGERQ2 dgerq2, DLARFB dlarfb, DLARFT dlarft, XERBLA xerbla, ILAENV ilaenv)
        {
    

            #region Set Dependencies
            
            _dgerq2 = dgerq2; _dlarfb = dlarfb; _dlarft = dlarft; _xerbla = xerbla; _ilaenv = ilaenv; 

            #endregion

        }
    
        public DGERQF()
        {
    

            #region Dependencies (Initialization)
            
            var lsame = new LSAME();
            var xerbla = new XERBLA();
            var dlamc3 = new DLAMC3();
            var dlapy2 = new DLAPY2();
            var dnrm2 = new DNRM2();
            var dscal = new DSCAL();
            var dcopy = new DCOPY();
            var ieeeck = new IEEECK();
            var iparmq = new IPARMQ();
            var dgemv = new DGEMV(lsame, xerbla);
            var dger = new DGER(xerbla);
            var dlarf = new DLARF(dgemv, dger, lsame);
            var dlamc1 = new DLAMC1(dlamc3);
            var dlamc4 = new DLAMC4(dlamc3);
            var dlamc5 = new DLAMC5(dlamc3);
            var dlamc2 = new DLAMC2(dlamc3, dlamc1, dlamc4, dlamc5);
            var dlamch = new DLAMCH(lsame, dlamc2);
            var dlarfg = new DLARFG(dlamch, dlapy2, dnrm2, dscal);
            var dgerq2 = new DGERQ2(dlarf, dlarfg, xerbla);
            var dgemm = new DGEMM(lsame, xerbla);
            var dtrmm = new DTRMM(lsame, xerbla);
            var dlarfb = new DLARFB(lsame, dcopy, dgemm, dtrmm);
            var dtrmv = new DTRMV(lsame, xerbla);
            var dlarft = new DLARFT(dgemv, dtrmv, lsame);
            var ilaenv = new ILAENV(ieeeck, iparmq);

            #endregion


            #region Set Dependencies
            
            _dgerq2 = dgerq2; _dlarfb = dlarfb; _dlarft = dlarft; _xerbla = xerbla; _ilaenv = ilaenv; 

            #endregion

        }
        /// <summary>
        /// Purpose
        /// =======
        /// 
        /// DGERQF computes an RQ factorization of a real M-by-N matrix A:
        /// A = R * Q.
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
        /// <param name="A">
        /// (input/output) DOUBLE PRECISION array, dimension (LDA,N)
        /// On entry, the M-by-N matrix A.
        /// On exit,
        /// if m .LE. n, the upper triangle of the subarray
        /// A(1:m,n-m+1:n) contains the M-by-M upper triangular matrix R;
        /// if m .GE. n, the elements on and above the (m-n)-th subdiagonal
        /// contain the M-by-N upper trapezoidal matrix R;
        /// the remaining elements, with the array TAU, represent the
        /// orthogonal matrix Q as a product of min(m,n) elementary
        /// reflectors (see Further Details).
        ///</param>
        /// <param name="LDA">
        /// (input) INTEGER
        /// The leading dimension of the array A.  LDA .GE. max(1,M).
        ///</param>
        /// <param name="TAU">
        /// (output) DOUBLE PRECISION array, dimension (min(M,N))
        /// The scalar factors of the elementary reflectors (see Further
        /// Details).
        ///</param>
        /// <param name="WORK">
        /// (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
        /// On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
        ///</param>
        /// <param name="LWORK">
        /// (input) INTEGER
        /// The dimension of the array WORK.  LWORK .GE. max(1,M).
        /// For optimum performance LWORK .GE. M*NB, where NB is
        /// the optimal blocksize.
        /// 
        /// If LWORK = -1, then a workspace query is assumed; the routine
        /// only calculates the optimal size of the WORK array, returns
        /// this value as the first entry of the WORK array, and no error
        /// message related to LWORK is issued by XERBLA.
        ///</param>
        /// <param name="INFO">
        /// (output) INTEGER
        /// = 0:  successful exit
        /// .LT. 0:  if INFO = -i, the i-th argument had an illegal value
        ///</param>
        public void Run(int M, int N, ref double[] A, int offset_a, int LDA, ref double[] TAU, int offset_tau, ref double[] WORK, int offset_work
                         , int LWORK, ref int INFO)
        {

            #region Variables
            
            var LQUERY = false; var I = 0; var IB = 0; var IINFO = 0; var IWS = 0; var K = 0; var KI = 0; var KK = 0; 
            var LDWORK = 0;var LWKOPT = 0; var MU = 0; var NB = 0; var NBMIN = 0; var NU = 0; var NX = 0; 

            #endregion


            #region Array Index Correction
            
             var o_a = -1 - LDA + offset_a;  var o_tau = -1 + offset_tau;  var o_work = -1 + offset_work; 

            #endregion


            #region Prolog
            
            // *
            // *  -- LAPACK routine (version 3.1) --
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
            // *  DGERQF computes an RQ factorization of a real M-by-N matrix A:
            // *  A = R * Q.
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
            // *  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
            // *          On entry, the M-by-N matrix A.
            // *          On exit,
            // *          if m <= n, the upper triangle of the subarray
            // *          A(1:m,n-m+1:n) contains the M-by-M upper triangular matrix R;
            // *          if m >= n, the elements on and above the (m-n)-th subdiagonal
            // *          contain the M-by-N upper trapezoidal matrix R;
            // *          the remaining elements, with the array TAU, represent the
            // *          orthogonal matrix Q as a product of min(m,n) elementary
            // *          reflectors (see Further Details).
            // *
            // *  LDA     (input) INTEGER
            // *          The leading dimension of the array A.  LDA >= max(1,M).
            // *
            // *  TAU     (output) DOUBLE PRECISION array, dimension (min(M,N))
            // *          The scalar factors of the elementary reflectors (see Further
            // *          Details).
            // *
            // *  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
            // *          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
            // *
            // *  LWORK   (input) INTEGER
            // *          The dimension of the array WORK.  LWORK >= max(1,M).
            // *          For optimum performance LWORK >= M*NB, where NB is
            // *          the optimal blocksize.
            // *
            // *          If LWORK = -1, then a workspace query is assumed; the routine
            // *          only calculates the optimal size of the WORK array, returns
            // *          this value as the first entry of the WORK array, and no error
            // *          message related to LWORK is issued by XERBLA.
            // *
            // *  INFO    (output) INTEGER
            // *          = 0:  successful exit
            // *          < 0:  if INFO = -i, the i-th argument had an illegal value
            // *
            // *  Further Details
            // *  ===============
            // *
            // *  The matrix Q is represented as a product of elementary reflectors
            // *
            // *     Q = H(1) H(2) . . . H(k), where k = min(m,n).
            // *
            // *  Each H(i) has the form
            // *
            // *     H(i) = I - tau * v * v'
            // *
            // *  where tau is a real scalar, and v is a real vector with
            // *  v(n-k+i+1:n) = 0 and v(n-k+i) = 1; v(1:n-k+i-1) is stored on exit in
            // *  A(m-k+i,1:n-k+i-1), and tau in TAU(i).
            // *
            // *  =====================================================================
            // *
            // *     .. Local Scalars ..
            // *     ..
            // *     .. External Subroutines ..
            // *     ..
            // *     .. Intrinsic Functions ..
            //      INTRINSIC          MAX, MIN;
            // *     ..
            // *     .. External Functions ..
            // *     ..
            // *     .. Executable Statements ..
            // *
            // *     Test the input arguments
            // *

            #endregion


            #region Body
            
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
                    if (LDA < Math.Max(1, M))
                    {
                        INFO =  - 4;
                    }
                }
            }
            // *
            if (INFO == 0)
            {
                K = Math.Min(M, N);
                if (K == 0)
                {
                    LWKOPT = 1;
                }
                else
                {
                    NB = _ilaenv.Run(1, "DGERQF", " ", M, N,  - 1,  - 1);
                    LWKOPT = M * NB;
                }
                WORK[1 + o_work] = LWKOPT;
                // *
                if (LWORK < Math.Max(1, M) && !LQUERY)
                {
                    INFO =  - 7;
                }
            }
            // *
            if (INFO != 0)
            {
                _xerbla.Run("DGERQF",  - INFO);
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
            if (K == 0)
            {
                return;
            }
            // *
            NBMIN = 2;
            NX = 1;
            IWS = M;
            if (NB > 1 && NB < K)
            {
                // *
                // *        Determine when to cross over from blocked to unblocked code.
                // *
                NX = Math.Max(0, _ilaenv.Run(3, "DGERQF", " ", M, N,  - 1,  - 1));
                if (NX < K)
                {
                    // *
                    // *           Determine if workspace is large enough for blocked code.
                    // *
                    LDWORK = M;
                    IWS = LDWORK * NB;
                    if (LWORK < IWS)
                    {
                        // *
                        // *              Not enough workspace to use optimal NB:  reduce NB and
                        // *              determine the minimum value of NB.
                        // *
                        NB = LWORK / LDWORK;
                        NBMIN = Math.Max(2, _ilaenv.Run(2, "DGERQF", " ", M, N,  - 1,  - 1));
                    }
                }
            }
            // *
            if (NB >= NBMIN && NB < K && NX < K)
            {
                // *
                // *        Use blocked code initially.
                // *        The last kk rows are handled by the block method.
                // *
                KI = (K - NX - 1) / NB * NB;
                KK = Math.Min(K, KI + NB);
                // *
                for (I = K - KK + KI + 1;  - NB >= 0 ? I <= K - KK + 1 : I >= K - KK + 1; I +=  - NB)
                {
                    IB = Math.Min(K - I + 1, NB);
                    // *
                    // *           Compute the RQ factorization of the current block
                    // *           A(m-k+i:m-k+i+ib-1,1:n-k+i+ib-1)
                    // *
                    _dgerq2.Run(IB, N - K + I + IB - 1, ref A, M - K + I+1 * LDA + o_a, LDA, ref TAU, I + o_tau, ref WORK, offset_work
                                     , ref IINFO);
                    if (M - K + I > 1)
                    {
                        // *
                        // *              Form the triangular factor of the block reflector
                        // *              H = H(i+ib-1) . . . H(i+1) H(i)
                        // *
                        _dlarft.Run("Backward", "Rowwise", N - K + I + IB - 1, IB, ref A, M - K + I+1 * LDA + o_a, LDA
                                         , TAU, I + o_tau, ref WORK, offset_work, LDWORK);
                        // *
                        // *              Apply H to A(1:m-k+i-1,1:n-k+i+ib-1) from the right
                        // *
                        _dlarfb.Run("Right", "No transpose", "Backward", "Rowwise", M - K + I - 1, N - K + I + IB - 1
                                         , IB, A, M - K + I+1 * LDA + o_a, LDA, WORK, offset_work, LDWORK, ref A, offset_a
                                         , LDA, ref WORK, IB + 1 + o_work, LDWORK);
                    }
                }
                MU = M - K + I + NB - 1;
                NU = N - K + I + NB - 1;
            }
            else
            {
                MU = M;
                NU = N;
            }
            // *
            // *     Use unblocked code to factor the last or only block
            // *
            if (MU > 0 && NU > 0)
            {
                _dgerq2.Run(MU, NU, ref A, offset_a, LDA, ref TAU, offset_tau, ref WORK, offset_work
                                 , ref IINFO);
            }
            // *
            WORK[1 + o_work] = IWS;
            return;
            // *
            // *     End of DGERQF
            // *

            #endregion

        }
    }
}
