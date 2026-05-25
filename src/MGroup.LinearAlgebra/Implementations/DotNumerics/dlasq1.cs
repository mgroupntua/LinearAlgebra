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
    /// DLASQ1 computes the singular values of a real N-by-N bidiagonal
    /// matrix with diagonal D and off-diagonal E. The singular values
    /// are computed to high relative accuracy, in the absence of
    /// denormalization, underflow and overflow. The algorithm was first
    /// presented in
    /// 
    /// "Accurate singular values and differential qd algorithms" by K. V.
    /// Fernando and B. N. Parlett, Numer. Math., Vol-67, No. 2, pp. 191-230,
    /// 1994,
    /// 
    /// and the present implementation is described in "An implementation of
    /// the dqds Algorithm (Positive Case)", LAPACK Working Note.
    /// 
    ///</summary>
    public class DLASQ1
    {
    

        #region Dependencies
        
        DCOPY _dcopy; DLAS2 _dlas2; DLASCL _dlascl; DLASQ2 _dlasq2; DLASRT _dlasrt; XERBLA _xerbla; DLAMCH _dlamch; 

        #endregion


        #region Variables
        
        const double ZERO = 0.0E0; 

        #endregion

        public DLASQ1(DCOPY dcopy, DLAS2 dlas2, DLASCL dlascl, DLASQ2 dlasq2, DLASRT dlasrt, XERBLA xerbla, DLAMCH dlamch)
        {
    

            #region Set Dependencies
            
            _dcopy = dcopy; _dlas2 = dlas2; _dlascl = dlascl; _dlasq2 = dlasq2; _dlasrt = dlasrt; 
            _xerbla = xerbla;_dlamch = dlamch; 

            #endregion

        }
    
        public DLASQ1()
        {
    

            #region Dependencies (Initialization)
            
            var dcopy = new DCOPY();
            var dlas2 = new DLAS2();
            var lsame = new LSAME();
            var dlamc3 = new DLAMC3();
            var xerbla = new XERBLA();
            var dlasq5 = new DLASQ5();
            var dlazq4 = new DLAZQ4();
            var ieeeck = new IEEECK();
            var iparmq = new IPARMQ();
            var dlamc1 = new DLAMC1(dlamc3);
            var dlamc4 = new DLAMC4(dlamc3);
            var dlamc5 = new DLAMC5(dlamc3);
            var dlamc2 = new DLAMC2(dlamc3, dlamc1, dlamc4, dlamc5);
            var dlamch = new DLAMCH(lsame, dlamc2);
            var dlascl = new DLASCL(lsame, dlamch, xerbla);
            var dlasq6 = new DLASQ6(dlamch);
            var dlazq3 = new DLAZQ3(dlasq5, dlasq6, dlazq4, dlamch);
            var dlasrt = new DLASRT(lsame, xerbla);
            var ilaenv = new ILAENV(ieeeck, iparmq);
            var dlasq2 = new DLASQ2(dlazq3, dlasrt, xerbla, dlamch, ilaenv);

            #endregion


            #region Set Dependencies
            
            _dcopy = dcopy; _dlas2 = dlas2; _dlascl = dlascl; _dlasq2 = dlasq2; _dlasrt = dlasrt; 
            _xerbla = xerbla;_dlamch = dlamch; 

            #endregion

        }
        /// <summary>
        /// Purpose
        /// =======
        /// 
        /// DLASQ1 computes the singular values of a real N-by-N bidiagonal
        /// matrix with diagonal D and off-diagonal E. The singular values
        /// are computed to high relative accuracy, in the absence of
        /// denormalization, underflow and overflow. The algorithm was first
        /// presented in
        /// 
        /// "Accurate singular values and differential qd algorithms" by K. V.
        /// Fernando and B. N. Parlett, Numer. Math., Vol-67, No. 2, pp. 191-230,
        /// 1994,
        /// 
        /// and the present implementation is described in "An implementation of
        /// the dqds Algorithm (Positive Case)", LAPACK Working Note.
        /// 
        ///</summary>
        /// <param name="N">
        /// (input) INTEGER
        /// The number of rows and columns in the matrix. N .GE. 0.
        ///</param>
        /// <param name="D">
        /// (input/output) DOUBLE PRECISION array, dimension (N)
        /// On entry, D contains the diagonal elements of the
        /// bidiagonal matrix whose SVD is desired. On normal exit,
        /// D contains the singular values in decreasing order.
        ///</param>
        /// <param name="E">
        /// (input/output) DOUBLE PRECISION array, dimension (N)
        /// On entry, elements E(1:N-1) contain the off-diagonal elements
        /// of the bidiagonal matrix whose SVD is desired.
        /// On exit, E is overwritten.
        ///</param>
        /// <param name="WORK">
        /// (workspace) DOUBLE PRECISION array, dimension (4*N)
        ///</param>
        /// <param name="INFO">
        /// (output) INTEGER
        /// = 0: successful exit
        /// .LT. 0: if INFO = -i, the i-th argument had an illegal value
        /// .GT. 0: the algorithm failed
        /// = 1, a split was marked by a positive value in E
        /// = 2, current block of Z not diagonalized after 30*N
        /// iterations (in inner while loop)
        /// = 3, termination criterion of outer while loop not met 
        /// (program created more than N unreduced blocks)
        ///</param>
        public void Run(int N, ref double[] D, int offset_d, double[] E, int offset_e, ref double[] WORK, int offset_work, ref int INFO)
        {

            #region Variables
            
            var I = 0; var IINFO = 0; double EPS = 0; double SCALE = 0; double SAFMIN = 0; double SIGMN = 0; double SIGMX = 0; 

            #endregion


            #region Array Index Correction
            
             var o_d = -1 + offset_d;  var o_e = -1 + offset_e;  var o_work = -1 + offset_work; 

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
            // *  DLASQ1 computes the singular values of a real N-by-N bidiagonal
            // *  matrix with diagonal D and off-diagonal E. The singular values
            // *  are computed to high relative accuracy, in the absence of
            // *  denormalization, underflow and overflow. The algorithm was first
            // *  presented in
            // *
            // *  "Accurate singular values and differential qd algorithms" by K. V.
            // *  Fernando and B. N. Parlett, Numer. Math., Vol-67, No. 2, pp. 191-230,
            // *  1994,
            // *
            // *  and the present implementation is described in "An implementation of
            // *  the dqds Algorithm (Positive Case)", LAPACK Working Note.
            // *
            // *  Arguments
            // *  =========
            // *
            // *  N     (input) INTEGER
            // *        The number of rows and columns in the matrix. N >= 0.
            // *
            // *  D     (input/output) DOUBLE PRECISION array, dimension (N)
            // *        On entry, D contains the diagonal elements of the
            // *        bidiagonal matrix whose SVD is desired. On normal exit,
            // *        D contains the singular values in decreasing order.
            // *
            // *  E     (input/output) DOUBLE PRECISION array, dimension (N)
            // *        On entry, elements E(1:N-1) contain the off-diagonal elements
            // *        of the bidiagonal matrix whose SVD is desired.
            // *        On exit, E is overwritten.
            // *
            // *  WORK  (workspace) DOUBLE PRECISION array, dimension (4*N)
            // *
            // *  INFO  (output) INTEGER
            // *        = 0: successful exit
            // *        < 0: if INFO = -i, the i-th argument had an illegal value
            // *        > 0: the algorithm failed
            // *             = 1, a split was marked by a positive value in E
            // *             = 2, current block of Z not diagonalized after 30*N
            // *                  iterations (in inner while loop)
            // *             = 3, termination criterion of outer while loop not met 
            // *                  (program created more than N unreduced blocks)
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
            //      INTRINSIC          ABS, MAX, SQRT;
            // *     ..
            // *     .. Executable Statements ..
            // *

            #endregion


            #region Body
            
            INFO = 0;
            if (N < 0)
            {
                INFO =  - 2;
                _xerbla.Run("DLASQ1",  - INFO);
                return;
            }
            else
            {
                if (N == 0)
                {
                    return;
                }
                else
                {
                    if (N == 1)
                    {
                        D[1 + o_d] = Math.Abs(D[1 + o_d]);
                        return;
                    }
                    else
                    {
                        if (N == 2)
                        {
                            _dlas2.Run(D[1 + o_d], E[1 + o_e], D[2 + o_d], ref SIGMN, ref SIGMX);
                            D[1 + o_d] = SIGMX;
                            D[2 + o_d] = SIGMN;
                            return;
                        }
                    }
                }
            }
            // *
            // *     Estimate the largest singular value.
            // *
            SIGMX = ZERO;
            for (I = 1; I <= N - 1; I++)
            {
                D[I + o_d] = Math.Abs(D[I + o_d]);
                SIGMX = Math.Max(SIGMX, Math.Abs(E[I + o_e]));
            }
            D[N + o_d] = Math.Abs(D[N + o_d]);
            // *
            // *     Early return if SIGMX is zero (matrix is already diagonal).
            // *
            if (SIGMX == ZERO)
            {
                _dlasrt.Run("D", N, ref D, offset_d, ref IINFO);
                return;
            }
            // *
            for (I = 1; I <= N; I++)
            {
                SIGMX = Math.Max(SIGMX, D[I + o_d]);
            }
            // *
            // *     Copy D and E into WORK (in the Z format) and scale (squaring the
            // *     input data makes scaling by a power of the radix pointless).
            // *
            EPS = _dlamch.Run("Precision");
            SAFMIN = _dlamch.Run("Safe minimum");
            SCALE = Math.Sqrt(EPS / SAFMIN);
            _dcopy.Run(N, D, offset_d, 1, ref WORK, 1 + o_work, 2);
            _dcopy.Run(N - 1, E, offset_e, 1, ref WORK, 2 + o_work, 2);
            _dlascl.Run("G", 0, 0, SIGMX, SCALE, 2 * N - 1
                             , 1, ref WORK, offset_work, 2 * N - 1, ref IINFO);
            // *         
            // *     Compute the q's and e's.
            // *
            for (I = 1; I <= 2 * N - 1; I++)
            {
                WORK[I + o_work] = Math.Pow(WORK[I + o_work],2);
            }
            WORK[2 * N + o_work] = ZERO;
            // *
            _dlasq2.Run(N, ref WORK, offset_work, ref INFO);
            // *
            if (INFO == 0)
            {
                for (I = 1; I <= N; I++)
                {
                    D[I + o_d] = Math.Sqrt(WORK[I + o_work]);
                }
                _dlascl.Run("G", 0, 0, SCALE, SIGMX, N
                                 , 1, ref D, offset_d, N, ref IINFO);
            }
            // *
            return;
            // *
            // *     End of DLASQ1
            // *

            #endregion

        }
    }
}
