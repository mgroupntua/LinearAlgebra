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
    /// DORGHR generates a real orthogonal matrix Q which is defined as the
    /// product of IHI-ILO elementary reflectors of order N, as returned by
    /// DGEHRD:
    /// 
    /// Q = H(ilo) H(ilo+1) . . . H(ihi-1).
    /// 
    ///</summary>
    public class DORGHR
    {
    

        #region Dependencies
        
        DORGQR _dorgqr; XERBLA _xerbla; ILAENV _ilaenv; 

        #endregion


        #region Variables
        
        const double ZERO = 0.0E+0; const double ONE = 1.0E+0; 

        #endregion

        public DORGHR(DORGQR dorgqr, XERBLA xerbla, ILAENV ilaenv)
        {
    

            #region Set Dependencies
            
            _dorgqr = dorgqr; _xerbla = xerbla; _ilaenv = ilaenv; 

            #endregion

        }
    
        public DORGHR()
        {
    

            #region Dependencies (Initialization)
            
            var lsame = new LSAME();
            var dcopy = new DCOPY();
            var xerbla = new XERBLA();
            var dscal = new DSCAL();
            var ieeeck = new IEEECK();
            var iparmq = new IPARMQ();
            var dgemm = new DGEMM(lsame, xerbla);
            var dtrmm = new DTRMM(lsame, xerbla);
            var dlarfb = new DLARFB(lsame, dcopy, dgemm, dtrmm);
            var dgemv = new DGEMV(lsame, xerbla);
            var dtrmv = new DTRMV(lsame, xerbla);
            var dlarft = new DLARFT(dgemv, dtrmv, lsame);
            var dger = new DGER(xerbla);
            var dlarf = new DLARF(dgemv, dger, lsame);
            var dorg2r = new DORG2R(dlarf, dscal, xerbla);
            var ilaenv = new ILAENV(ieeeck, iparmq);
            var dorgqr = new DORGQR(dlarfb, dlarft, dorg2r, xerbla, ilaenv);

            #endregion


            #region Set Dependencies
            
            _dorgqr = dorgqr; _xerbla = xerbla; _ilaenv = ilaenv; 

            #endregion

        }
        /// <summary>
        /// Purpose
        /// =======
        /// 
        /// DORGHR generates a real orthogonal matrix Q which is defined as the
        /// product of IHI-ILO elementary reflectors of order N, as returned by
        /// DGEHRD:
        /// 
        /// Q = H(ilo) H(ilo+1) . . . H(ihi-1).
        /// 
        ///</summary>
        /// <param name="N">
        /// (input) INTEGER
        /// The order of the matrix Q. N .GE. 0.
        ///</param>
        /// <param name="ILO">
        /// (input) INTEGER
        ///</param>
        /// <param name="IHI">
        /// (input) INTEGER
        /// ILO and IHI must have the same values as in the previous call
        /// of DGEHRD. Q is equal to the unit matrix except in the
        /// submatrix Q(ilo+1:ihi,ilo+1:ihi).
        /// 1 .LE. ILO .LE. IHI .LE. N, if N .GT. 0; ILO=1 and IHI=0, if N=0.
        ///</param>
        /// <param name="A">
        /// (input/output) DOUBLE PRECISION array, dimension (LDA,N)
        /// On entry, the vectors which define the elementary reflectors,
        /// as returned by DGEHRD.
        /// On exit, the N-by-N orthogonal matrix Q.
        ///</param>
        /// <param name="LDA">
        /// (input) INTEGER
        /// The leading dimension of the array A. LDA .GE. max(1,N).
        ///</param>
        /// <param name="TAU">
        /// (input) DOUBLE PRECISION array, dimension (N-1)
        /// TAU(i) must contain the scalar factor of the elementary
        /// reflector H(i), as returned by DGEHRD.
        ///</param>
        /// <param name="WORK">
        /// (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
        /// On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
        ///</param>
        /// <param name="LWORK">
        /// (input) INTEGER
        /// The dimension of the array WORK. LWORK .GE. IHI-ILO.
        /// For optimum performance LWORK .GE. (IHI-ILO)*NB, where NB is
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
        public void Run(int N, int ILO, int IHI, ref double[] A, int offset_a, int LDA, double[] TAU, int offset_tau
                         , ref double[] WORK, int offset_work, int LWORK, ref int INFO)
        {

            #region Variables
            
            var LQUERY = false; var I = 0; var IINFO = 0; var J = 0; var LWKOPT = 0; var NB = 0; var NH = 0; 

            #endregion


            #region Implicit Variables
            
            var A_J = 0; var A_0 = 0; 

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
            // *  DORGHR generates a real orthogonal matrix Q which is defined as the
            // *  product of IHI-ILO elementary reflectors of order N, as returned by
            // *  DGEHRD:
            // *
            // *  Q = H(ilo) H(ilo+1) . . . H(ihi-1).
            // *
            // *  Arguments
            // *  =========
            // *
            // *  N       (input) INTEGER
            // *          The order of the matrix Q. N >= 0.
            // *
            // *  ILO     (input) INTEGER
            // *  IHI     (input) INTEGER
            // *          ILO and IHI must have the same values as in the previous call
            // *          of DGEHRD. Q is equal to the unit matrix except in the
            // *          submatrix Q(ilo+1:ihi,ilo+1:ihi).
            // *          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.
            // *
            // *  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
            // *          On entry, the vectors which define the elementary reflectors,
            // *          as returned by DGEHRD.
            // *          On exit, the N-by-N orthogonal matrix Q.
            // *
            // *  LDA     (input) INTEGER
            // *          The leading dimension of the array A. LDA >= max(1,N).
            // *
            // *  TAU     (input) DOUBLE PRECISION array, dimension (N-1)
            // *          TAU(i) must contain the scalar factor of the elementary
            // *          reflector H(i), as returned by DGEHRD.
            // *
            // *  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
            // *          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
            // *
            // *  LWORK   (input) INTEGER
            // *          The dimension of the array WORK. LWORK >= IHI-ILO.
            // *          For optimum performance LWORK >= (IHI-ILO)*NB, where NB is
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
            // *  =====================================================================
            // *
            // *     .. Parameters ..
            // *     ..
            // *     .. Local Scalars ..
            // *     ..
            // *     .. External Subroutines ..
            // *     ..
            // *     .. External Functions ..
            // *     ..
            // *     .. Intrinsic Functions ..
            //      INTRINSIC          MAX, MIN;
            // *     ..
            // *     .. Executable Statements ..
            // *
            // *     Test the input arguments
            // *

            #endregion


            #region Body
            
            INFO = 0;
            NH = IHI - ILO;
            LQUERY = LWORK ==  - 1;
            if (N < 0)
            {
                INFO =  - 1;
            }
            else
            {
                if (ILO < 1 || ILO > Math.Max(1, N))
                {
                    INFO =  - 2;
                }
                else
                {
                    if (IHI < Math.Min(ILO, N) || IHI > N)
                    {
                        INFO =  - 3;
                    }
                    else
                    {
                        if (LDA < Math.Max(1, N))
                        {
                            INFO =  - 5;
                        }
                        else
                        {
                            if (LWORK < Math.Max(1, NH) && !LQUERY)
                            {
                                INFO =  - 8;
                            }
                        }
                    }
                }
            }
            // *
            if (INFO == 0)
            {
                NB = _ilaenv.Run(1, "DORGQR", " ", NH, NH, NH,  - 1);
                LWKOPT = Math.Max(1, NH) * NB;
                WORK[1 + o_work] = LWKOPT;
            }
            // *
            if (INFO != 0)
            {
                _xerbla.Run("DORGHR",  - INFO);
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
            if (N == 0)
            {
                WORK[1 + o_work] = 1;
                return;
            }
            // *
            // *     Shift the vectors which define the elementary reflectors one
            // *     column to the right, and set the first ilo and the last n-ihi
            // *     rows and columns to those of the unit matrix
            // *
            for (J = IHI; J >= ILO + 1; J +=  - 1)
            {
                A_J = J * LDA + o_a;
                for (I = 1; I <= J - 1; I++)
                {
                    A[I + A_J] = ZERO;
                }
                A_J = J * LDA + o_a;
                A_0 = (J - 1) * LDA + o_a;
                for (I = J + 1; I <= IHI; I++)
                {
                    A[I + A_J] = A[I + A_0];
                }
                A_J = J * LDA + o_a;
                for (I = IHI + 1; I <= N; I++)
                {
                    A[I + A_J] = ZERO;
                }
            }
            for (J = 1; J <= ILO; J++)
            {
                A_J = J * LDA + o_a;
                for (I = 1; I <= N; I++)
                {
                    A[I + A_J] = ZERO;
                }
                A[J+J * LDA + o_a] = ONE;
            }
            for (J = IHI + 1; J <= N; J++)
            {
                A_J = J * LDA + o_a;
                for (I = 1; I <= N; I++)
                {
                    A[I + A_J] = ZERO;
                }
                A[J+J * LDA + o_a] = ONE;
            }
            // *
            if (NH > 0)
            {
                // *
                // *        Generate Q(ilo+1:ihi,ilo+1:ihi)
                // *
                _dorgqr.Run(NH, NH, NH, ref A, ILO + 1+(ILO + 1) * LDA + o_a, LDA, TAU, ILO + o_tau
                                 , ref WORK, offset_work, LWORK, ref IINFO);
            }
            WORK[1 + o_work] = LWKOPT;
            return;
            // *
            // *     End of DORGHR
            // *

            #endregion

        }
    }
}
