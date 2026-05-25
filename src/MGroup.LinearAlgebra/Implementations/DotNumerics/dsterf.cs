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

using MGroup.LinearAlgebra.Implementations.DotNumerics.FortranLibrary;

namespace MGroup.LinearAlgebra.Implementations.DotNumerics
{
    /// <summary>
    /// -- LAPACK routine (version 3.1) --
    /// Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
    /// November 2006
    /// Purpose
    /// =======
    /// 
    /// DSTERF computes all eigenvalues of a symmetric tridiagonal matrix
    /// using the Pal-Walker-Kahan variant of the QL or QR algorithm.
    /// 
    ///</summary>
    public class DSTERF
    {
    

        #region Dependencies
        
        DLAMCH _dlamch; DLANST _dlanst; DLAPY2 _dlapy2; DLAE2 _dlae2; DLASCL _dlascl; DLASRT _dlasrt; XERBLA _xerbla; 

        #endregion


        #region Variables
        
        const double ZERO = 0.0E0; const double ONE = 1.0E0; const double TWO = 2.0E0; const double THREE = 3.0E0; 
        const int MAXIT = 30;

        #endregion

        public DSTERF(DLAMCH dlamch, DLANST dlanst, DLAPY2 dlapy2, DLAE2 dlae2, DLASCL dlascl, DLASRT dlasrt, XERBLA xerbla)
        {
    

            #region Set Dependencies
            
            _dlamch = dlamch; _dlanst = dlanst; _dlapy2 = dlapy2; _dlae2 = dlae2; _dlascl = dlascl; 
            _dlasrt = dlasrt;_xerbla = xerbla; 

            #endregion

        }
    
        public DSTERF()
        {
    

            #region Dependencies (Initialization)
            
            var lsame = new LSAME();
            var dlamc3 = new DLAMC3();
            var dlassq = new DLASSQ();
            var dlapy2 = new DLAPY2();
            var dlae2 = new DLAE2();
            var xerbla = new XERBLA();
            var dlamc1 = new DLAMC1(dlamc3);
            var dlamc4 = new DLAMC4(dlamc3);
            var dlamc5 = new DLAMC5(dlamc3);
            var dlamc2 = new DLAMC2(dlamc3, dlamc1, dlamc4, dlamc5);
            var dlamch = new DLAMCH(lsame, dlamc2);
            var dlanst = new DLANST(lsame, dlassq);
            var dlascl = new DLASCL(lsame, dlamch, xerbla);
            var dlasrt = new DLASRT(lsame, xerbla);

            #endregion


            #region Set Dependencies
            
            _dlamch = dlamch; _dlanst = dlanst; _dlapy2 = dlapy2; _dlae2 = dlae2; _dlascl = dlascl; 
            _dlasrt = dlasrt;_xerbla = xerbla; 

            #endregion

        }
        /// <summary>
        /// Purpose
        /// =======
        /// 
        /// DSTERF computes all eigenvalues of a symmetric tridiagonal matrix
        /// using the Pal-Walker-Kahan variant of the QL or QR algorithm.
        /// 
        ///</summary>
        /// <param name="N">
        /// (input) INTEGER
        /// The order of the matrix.  N .GE. 0.
        ///</param>
        /// <param name="D">
        /// (input/output) DOUBLE PRECISION array, dimension (N)
        /// On entry, the n diagonal elements of the tridiagonal matrix.
        /// On exit, if INFO = 0, the eigenvalues in ascending order.
        ///</param>
        /// <param name="E">
        /// (input/output) DOUBLE PRECISION array, dimension (N-1)
        /// On entry, the (n-1) subdiagonal elements of the tridiagonal
        /// matrix.
        /// On exit, E has been destroyed.
        ///</param>
        /// <param name="INFO">
        /// (output) INTEGER
        /// = 0:  successful exit
        /// .LT. 0:  if INFO = -i, the i-th argument had an illegal value
        /// .GT. 0:  the algorithm failed to find all of the eigenvalues in
        /// a total of 30*N iterations; if INFO = i, then i
        /// elements of E have not converged to zero.
        ///</param>
        public void Run(int N, ref double[] D, int offset_d, ref double[] E, int offset_e, ref int INFO)
        {

            #region Variables
            
            var I = 0; var ISCALE = 0; var JTOT = 0; var L = 0; var L1 = 0; var LEND = 0; var LENDSV = 0; var LSV = 0; var M = 0; 
            var NMAXIT = 0;double ALPHA = 0; double ANORM = 0; double BB = 0; double C = 0; double EPS = 0; double EPS2 = 0; 
            double GAMMA = 0;double OLDC = 0; double OLDGAM = 0; double P = 0; double R = 0; double RT1 = 0; double RT2 = 0; 
            double RTE = 0;double S = 0; double SAFMAX = 0; double SAFMIN = 0; double SIGMA = 0; double SSFMAX = 0; 
            double SSFMIN = 0;

            #endregion


            #region Array Index Correction
            
             var o_d = -1 + offset_d;  var o_e = -1 + offset_e; 

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
            // *  DSTERF computes all eigenvalues of a symmetric tridiagonal matrix
            // *  using the Pal-Walker-Kahan variant of the QL or QR algorithm.
            // *
            // *  Arguments
            // *  =========
            // *
            // *  N       (input) INTEGER
            // *          The order of the matrix.  N >= 0.
            // *
            // *  D       (input/output) DOUBLE PRECISION array, dimension (N)
            // *          On entry, the n diagonal elements of the tridiagonal matrix.
            // *          On exit, if INFO = 0, the eigenvalues in ascending order.
            // *
            // *  E       (input/output) DOUBLE PRECISION array, dimension (N-1)
            // *          On entry, the (n-1) subdiagonal elements of the tridiagonal
            // *          matrix.
            // *          On exit, E has been destroyed.
            // *
            // *  INFO    (output) INTEGER
            // *          = 0:  successful exit
            // *          < 0:  if INFO = -i, the i-th argument had an illegal value
            // *          > 0:  the algorithm failed to find all of the eigenvalues in
            // *                a total of 30*N iterations; if INFO = i, then i
            // *                elements of E have not converged to zero.
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
            //      INTRINSIC          ABS, SIGN, SQRT;
            // *     ..
            // *     .. Executable Statements ..
            // *
            // *     Test the input parameters.
            // *

            #endregion


            #region Body
            
            INFO = 0;
            // *
            // *     Quick return if possible
            // *
            if (N < 0)
            {
                INFO =  - 1;
                _xerbla.Run("DSTERF",  - INFO);
                return;
            }
            if (N <= 1) return;
            // *
            // *     Determine the unit roundoff for this environment.
            // *
            EPS = _dlamch.Run("E");
            EPS2 = Math.Pow(EPS,2);
            SAFMIN = _dlamch.Run("S");
            SAFMAX = ONE / SAFMIN;
            SSFMAX = Math.Sqrt(SAFMAX) / THREE;
            SSFMIN = Math.Sqrt(SAFMIN) / EPS2;
            // *
            // *     Compute the eigenvalues of the tridiagonal matrix.
            // *
            NMAXIT = N * MAXIT;
            SIGMA = ZERO;
            JTOT = 0;
            // *
            // *     Determine where the matrix splits and choose QL or QR iteration
            // *     for each block, according to whether top or bottom diagonal
            // *     element is smaller.
            // *
            L1 = 1;
            // *
        LABEL10:;
            if (L1 > N) goto LABEL170;
            if (L1 > 1) E[L1 - 1 + o_e] = ZERO;
            for (M = L1; M <= N - 1; M++)
            {
                if (Math.Abs(E[M + o_e]) <= Math.Sqrt(Math.Abs(D[M + o_d])) * Math.Sqrt(Math.Abs(D[M + 1 + o_d])) * EPS)
                {
                    E[M + o_e] = ZERO;
                    goto LABEL30;
                }
            }
            M = N;
            // *
        LABEL30:;
            L = L1;
            LSV = L;
            LEND = M;
            LENDSV = LEND;
            L1 = M + 1;
            if (LEND == L) goto LABEL10;
            // *
            // *     Scale submatrix in rows and columns L to LEND
            // *
            ANORM = _dlanst.Run("I", LEND - L + 1, D, L + o_d, E, L + o_e);
            ISCALE = 0;
            if (ANORM > SSFMAX)
            {
                ISCALE = 1;
                _dlascl.Run("G", 0, 0, ANORM, SSFMAX, LEND - L + 1
                                 , 1, ref D, L + o_d, N, ref INFO);
                _dlascl.Run("G", 0, 0, ANORM, SSFMAX, LEND - L
                                 , 1, ref E, L + o_e, N, ref INFO);
            }
            else
            {
                if (ANORM < SSFMIN)
                {
                    ISCALE = 2;
                    _dlascl.Run("G", 0, 0, ANORM, SSFMIN, LEND - L + 1
                                     , 1, ref D, L + o_d, N, ref INFO);
                    _dlascl.Run("G", 0, 0, ANORM, SSFMIN, LEND - L
                                     , 1, ref E, L + o_e, N, ref INFO);
                }
            }
            // *
            for (I = L; I <= LEND - 1; I++)
            {
                E[I + o_e] = Math.Pow(E[I + o_e],2);
            }
            // *
            // *     Choose between QL and QR iteration
            // *
            if (Math.Abs(D[LEND + o_d]) < Math.Abs(D[L + o_d]))
            {
                LEND = LSV;
                L = LENDSV;
            }
            // *
            if (LEND >= L)
            {
                // *
                // *        QL Iteration
                // *
                // *        Look for small subdiagonal element.
                // *
            LABEL50:;
                if (L != LEND)
                {
                    for (M = L; M <= LEND - 1; M++)
                    {
                        if (Math.Abs(E[M + o_e]) <= EPS2 * Math.Abs(D[M + o_d] * D[M + 1 + o_d])) goto LABEL70;
                    }
                }
                M = LEND;
                // *
            LABEL70:;
                if (M < LEND) E[M + o_e] = ZERO;
                P = D[L + o_d];
                if (M == L) goto LABEL90;
                // *
                // *        If remaining matrix is 2 by 2, use DLAE2 to compute its
                // *        eigenvalues.
                // *
                if (M == L + 1)
                {
                    RTE = Math.Sqrt(E[L + o_e]);
                    _dlae2.Run(D[L + o_d], RTE, D[L + 1 + o_d], ref RT1, ref RT2);
                    D[L + o_d] = RT1;
                    D[L + 1 + o_d] = RT2;
                    E[L + o_e] = ZERO;
                    L += 2;
                    if (L <= LEND) goto LABEL50;
                    goto LABEL150;
                }
                // *
                if (JTOT == NMAXIT) goto LABEL150;
                JTOT += 1;
                // *
                // *        Form shift.
                // *
                RTE = Math.Sqrt(E[L + o_e]);
                SIGMA = (D[L + 1 + o_d] - P) / (TWO * RTE);
                R = _dlapy2.Run(SIGMA, ONE);
                SIGMA = P - RTE / (SIGMA + FortranLib.Sign(R,SIGMA));
                // *
                C = ONE;
                S = ZERO;
                GAMMA = D[M + o_d] - SIGMA;
                P = GAMMA * GAMMA;
                // *
                // *        Inner loop
                // *
                for (I = M - 1; I >= L; I +=  - 1)
                {
                    BB = E[I + o_e];
                    R = P + BB;
                    if (I != M - 1) E[I + 1 + o_e] = S * R;
                    OLDC = C;
                    C = P / R;
                    S = BB / R;
                    OLDGAM = GAMMA;
                    ALPHA = D[I + o_d];
                    GAMMA = C * (ALPHA - SIGMA) - S * OLDGAM;
                    D[I + 1 + o_d] = OLDGAM + (ALPHA - GAMMA);
                    if (C != ZERO)
                    {
                        P = GAMMA * GAMMA / C;
                    }
                    else
                    {
                        P = OLDC * BB;
                    }
                }
                // *
                E[L + o_e] = S * P;
                D[L + o_d] = SIGMA + GAMMA;
                goto LABEL50;
                // *
                // *        Eigenvalue found.
                // *
            LABEL90:;
                D[L + o_d] = P;
                // *
                L += 1;
                if (L <= LEND) goto LABEL50;
                goto LABEL150;
                // *
            }
            else
            {
                // *
                // *        QR Iteration
                // *
                // *        Look for small superdiagonal element.
                // *
            LABEL100:;
                for (M = L; M >= LEND + 1; M +=  - 1)
                {
                    if (Math.Abs(E[M - 1 + o_e]) <= EPS2 * Math.Abs(D[M + o_d] * D[M - 1 + o_d])) goto LABEL120;
                }
                M = LEND;
                // *
            LABEL120:;
                if (M > LEND) E[M - 1 + o_e] = ZERO;
                P = D[L + o_d];
                if (M == L) goto LABEL140;
                // *
                // *        If remaining matrix is 2 by 2, use DLAE2 to compute its
                // *        eigenvalues.
                // *
                if (M == L - 1)
                {
                    RTE = Math.Sqrt(E[L - 1 + o_e]);
                    _dlae2.Run(D[L + o_d], RTE, D[L - 1 + o_d], ref RT1, ref RT2);
                    D[L + o_d] = RT1;
                    D[L - 1 + o_d] = RT2;
                    E[L - 1 + o_e] = ZERO;
                    L -= 2;
                    if (L >= LEND) goto LABEL100;
                    goto LABEL150;
                }
                // *
                if (JTOT == NMAXIT) goto LABEL150;
                JTOT += 1;
                // *
                // *        Form shift.
                // *
                RTE = Math.Sqrt(E[L - 1 + o_e]);
                SIGMA = (D[L - 1 + o_d] - P) / (TWO * RTE);
                R = _dlapy2.Run(SIGMA, ONE);
                SIGMA = P - RTE / (SIGMA + FortranLib.Sign(R,SIGMA));
                // *
                C = ONE;
                S = ZERO;
                GAMMA = D[M + o_d] - SIGMA;
                P = GAMMA * GAMMA;
                // *
                // *        Inner loop
                // *
                for (I = M; I <= L - 1; I++)
                {
                    BB = E[I + o_e];
                    R = P + BB;
                    if (I != M) E[I - 1 + o_e] = S * R;
                    OLDC = C;
                    C = P / R;
                    S = BB / R;
                    OLDGAM = GAMMA;
                    ALPHA = D[I + 1 + o_d];
                    GAMMA = C * (ALPHA - SIGMA) - S * OLDGAM;
                    D[I + o_d] = OLDGAM + (ALPHA - GAMMA);
                    if (C != ZERO)
                    {
                        P = GAMMA * GAMMA / C;
                    }
                    else
                    {
                        P = OLDC * BB;
                    }
                }
                // *
                E[L - 1 + o_e] = S * P;
                D[L + o_d] = SIGMA + GAMMA;
                goto LABEL100;
                // *
                // *        Eigenvalue found.
                // *
            LABEL140:;
                D[L + o_d] = P;
                // *
                L -= 1;
                if (L >= LEND) goto LABEL100;
                goto LABEL150;
                // *
            }
            // *
            // *     Undo scaling if necessary
            // *
        LABEL150:;
            if (ISCALE == 1)
            {
                _dlascl.Run("G", 0, 0, SSFMAX, ANORM, LENDSV - LSV + 1
                                 , 1, ref D, LSV + o_d, N, ref INFO);
            }
            if (ISCALE == 2)
            {
                _dlascl.Run("G", 0, 0, SSFMIN, ANORM, LENDSV - LSV + 1
                                 , 1, ref D, LSV + o_d, N, ref INFO);
            }
            // *
            // *     Check for no convergence to an eigenvalue after a total
            // *     of N*MAXIT iterations.
            // *
            if (JTOT < NMAXIT) goto LABEL10;
            for (I = 1; I <= N - 1; I++)
            {
                if (E[I + o_e] != ZERO) INFO += 1;
            }
            goto LABEL180;
            // *
            // *     Sort eigenvalues in increasing order.
            // *
        LABEL170:;
            _dlasrt.Run("I", N, ref D, offset_d, ref INFO);
            // *
        LABEL180:;
            return;
            // *
            // *     End of DSTERF
            // *

            #endregion

        }
    }
}
