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
    /// DSBEV computes all the eigenvalues and, optionally, eigenvectors of
    /// a real symmetric band matrix A.
    /// 
    ///</summary>
    public class DSBEV
    {
    

        #region Dependencies
        
        LSAME _lsame; DLAMCH _dlamch; DLANSB _dlansb; DLASCL _dlascl; DSBTRD _dsbtrd; DSCAL _dscal; DSTEQR _dsteqr; 
        DSTERF _dsterf;XERBLA _xerbla; 

        #endregion


        #region Variables
        
        const double ZERO = 0.0E0; const double ONE = 1.0E0; 

        #endregion

        public DSBEV(LSAME lsame, DLAMCH dlamch, DLANSB dlansb, DLASCL dlascl, DSBTRD dsbtrd, DSCAL dscal, DSTEQR dsteqr, DSTERF dsterf, XERBLA xerbla)
        {
    

            #region Set Dependencies
            
            _lsame = lsame; _dlamch = dlamch; _dlansb = dlansb; _dlascl = dlascl; _dsbtrd = dsbtrd; 
            _dscal = dscal;_dsteqr = dsteqr; _dsterf = dsterf; _xerbla = xerbla; 

            #endregion

        }
    
        public DSBEV()
        {
    

            #region Dependencies (Initialization)
            
            var lsame = new LSAME();
            var dlamc3 = new DLAMC3();
            var dlassq = new DLASSQ();
            var xerbla = new XERBLA();
            var dlar2v = new DLAR2V();
            var dlargv = new DLARGV();
            var dlartv = new DLARTV();
            var drot = new DROT();
            var dscal = new DSCAL();
            var dlapy2 = new DLAPY2();
            var dlae2 = new DLAE2();
            var dlaev2 = new DLAEV2();
            var dswap = new DSWAP();
            var dlamc1 = new DLAMC1(dlamc3);
            var dlamc4 = new DLAMC4(dlamc3);
            var dlamc5 = new DLAMC5(dlamc3);
            var dlamc2 = new DLAMC2(dlamc3, dlamc1, dlamc4, dlamc5);
            var dlamch = new DLAMCH(lsame, dlamc2);
            var dlansb = new DLANSB(dlassq, lsame);
            var dlascl = new DLASCL(lsame, dlamch, xerbla);
            var dlartg = new DLARTG(dlamch);
            var dlaset = new DLASET(lsame);
            var dsbtrd = new DSBTRD(dlar2v, dlargv, dlartg, dlartv, dlaset, drot, xerbla, lsame);
            var dlanst = new DLANST(lsame, dlassq);
            var dlasr = new DLASR(lsame, xerbla);
            var dlasrt = new DLASRT(lsame, xerbla);
            var dsteqr = new DSTEQR(lsame, dlamch, dlanst, dlapy2, dlae2, dlaev2, dlartg, dlascl, dlaset, dlasr
                                       , dlasrt, dswap, xerbla);
            var dsterf = new DSTERF(dlamch, dlanst, dlapy2, dlae2, dlascl, dlasrt, xerbla);

            #endregion


            #region Set Dependencies
            
            _lsame = lsame; _dlamch = dlamch; _dlansb = dlansb; _dlascl = dlascl; _dsbtrd = dsbtrd; 
            _dscal = dscal;_dsteqr = dsteqr; _dsterf = dsterf; _xerbla = xerbla; 

            #endregion

        }
        /// <summary>
        /// Purpose
        /// =======
        /// 
        /// DSBEV computes all the eigenvalues and, optionally, eigenvectors of
        /// a real symmetric band matrix A.
        /// 
        ///</summary>
        /// <param name="JOBZ">
        /// (input) CHARACTER*1
        /// = 'N':  Compute eigenvalues only;
        /// = 'V':  Compute eigenvalues and eigenvectors.
        ///</param>
        /// <param name="UPLO">
        /// (input) CHARACTER*1
        /// = 'U':  Upper triangle of A is stored;
        /// = 'L':  Lower triangle of A is stored.
        ///</param>
        /// <param name="N">
        /// (input) INTEGER
        /// The order of the matrix A.  N .GE. 0.
        ///</param>
        /// <param name="KD">
        /// (input) INTEGER
        /// The number of superdiagonals of the matrix A if UPLO = 'U',
        /// or the number of subdiagonals if UPLO = 'L'.  KD .GE. 0.
        ///</param>
        /// <param name="AB">
        /// (input/output) DOUBLE PRECISION array, dimension (LDAB, N)
        /// On entry, the upper or lower triangle of the symmetric band
        /// matrix A, stored in the first KD+1 rows of the array.  The
        /// j-th column of A is stored in the j-th column of the array AB
        /// as follows:
        /// if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd).LE.i.LE.j;
        /// if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j.LE.i.LE.min(n,j+kd).
        /// 
        /// On exit, AB is overwritten by values generated during the
        /// reduction to tridiagonal form.  If UPLO = 'U', the first
        /// superdiagonal and the diagonal of the tridiagonal matrix T
        /// are returned in rows KD and KD+1 of AB, and if UPLO = 'L',
        /// the diagonal and first subdiagonal of T are returned in the
        /// first two rows of AB.
        ///</param>
        /// <param name="LDAB">
        /// (input) INTEGER
        /// The leading dimension of the array AB.  LDAB .GE. KD + 1.
        ///</param>
        /// <param name="W">
        /// (output) DOUBLE PRECISION array, dimension (N)
        /// If INFO = 0, the eigenvalues in ascending order.
        ///</param>
        /// <param name="Z">
        /// (output) DOUBLE PRECISION array, dimension (LDZ, N)
        /// If JOBZ = 'V', then if INFO = 0, Z contains the orthonormal
        /// eigenvectors of the matrix A, with the i-th column of Z
        /// holding the eigenvector associated with W(i).
        /// If JOBZ = 'N', then Z is not referenced.
        ///</param>
        /// <param name="LDZ">
        /// (input) INTEGER
        /// The leading dimension of the array Z.  LDZ .GE. 1, and if
        /// JOBZ = 'V', LDZ .GE. max(1,N).
        ///</param>
        /// <param name="WORK">
        /// (workspace) DOUBLE PRECISION array, dimension (max(1,3*N-2))
        ///</param>
        /// <param name="INFO">
        /// (output) INTEGER
        /// = 0:  successful exit
        /// .LT. 0:  if INFO = -i, the i-th argument had an illegal value
        /// .GT. 0:  if INFO = i, the algorithm failed to converge; i
        /// off-diagonal elements of an intermediate tridiagonal
        /// form did not converge to zero.
        ///</param>
        public void Run(string JOBZ, string UPLO, int N, int KD, ref double[] AB, int offset_ab, int LDAB
                         , ref double[] W, int offset_w, ref double[] Z, int offset_z, int LDZ, ref double[] WORK, int offset_work, ref int INFO)
        {

            #region Variables
            
            var LOWER = false; var WANTZ = false; var IINFO = 0; var IMAX = 0; var INDE = 0; var INDWRK = 0; var ISCALE = 0; 
            double ANRM = 0;double BIGNUM = 0; double EPS = 0; double RMAX = 0; double RMIN = 0; double SAFMIN = 0; 
            double SIGMA = 0;double SMLNUM = 0; 

            #endregion


            #region Array Index Correction
            
             var o_ab = -1 - LDAB + offset_ab;  var o_w = -1 + offset_w;  var o_z = -1 - LDZ + offset_z; 
             var o_work = -1 + offset_work;

            #endregion


            #region Strings
            
            JOBZ = JOBZ.Substring(0, 1);  UPLO = UPLO.Substring(0, 1);  

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
            // *  DSBEV computes all the eigenvalues and, optionally, eigenvectors of
            // *  a real symmetric band matrix A.
            // *
            // *  Arguments
            // *  =========
            // *
            // *  JOBZ    (input) CHARACTER*1
            // *          = 'N':  Compute eigenvalues only;
            // *          = 'V':  Compute eigenvalues and eigenvectors.
            // *
            // *  UPLO    (input) CHARACTER*1
            // *          = 'U':  Upper triangle of A is stored;
            // *          = 'L':  Lower triangle of A is stored.
            // *
            // *  N       (input) INTEGER
            // *          The order of the matrix A.  N >= 0.
            // *
            // *  KD      (input) INTEGER
            // *          The number of superdiagonals of the matrix A if UPLO = 'U',
            // *          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.
            // *
            // *  AB      (input/output) DOUBLE PRECISION array, dimension (LDAB, N)
            // *          On entry, the upper or lower triangle of the symmetric band
            // *          matrix A, stored in the first KD+1 rows of the array.  The
            // *          j-th column of A is stored in the j-th column of the array AB
            // *          as follows:
            // *          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
            // *          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).
            // *
            // *          On exit, AB is overwritten by values generated during the
            // *          reduction to tridiagonal form.  If UPLO = 'U', the first
            // *          superdiagonal and the diagonal of the tridiagonal matrix T
            // *          are returned in rows KD and KD+1 of AB, and if UPLO = 'L',
            // *          the diagonal and first subdiagonal of T are returned in the
            // *          first two rows of AB.
            // *
            // *  LDAB    (input) INTEGER
            // *          The leading dimension of the array AB.  LDAB >= KD + 1.
            // *
            // *  W       (output) DOUBLE PRECISION array, dimension (N)
            // *          If INFO = 0, the eigenvalues in ascending order.
            // *
            // *  Z       (output) DOUBLE PRECISION array, dimension (LDZ, N)
            // *          If JOBZ = 'V', then if INFO = 0, Z contains the orthonormal
            // *          eigenvectors of the matrix A, with the i-th column of Z
            // *          holding the eigenvector associated with W(i).
            // *          If JOBZ = 'N', then Z is not referenced.
            // *
            // *  LDZ     (input) INTEGER
            // *          The leading dimension of the array Z.  LDZ >= 1, and if
            // *          JOBZ = 'V', LDZ >= max(1,N).
            // *
            // *  WORK    (workspace) DOUBLE PRECISION array, dimension (max(1,3*N-2))
            // *
            // *  INFO    (output) INTEGER
            // *          = 0:  successful exit
            // *          < 0:  if INFO = -i, the i-th argument had an illegal value
            // *          > 0:  if INFO = i, the algorithm failed to converge; i
            // *                off-diagonal elements of an intermediate tridiagonal
            // *                form did not converge to zero.
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
            //      INTRINSIC          SQRT;
            // *     ..
            // *     .. Executable Statements ..
            // *
            // *     Test the input parameters.
            // *

            #endregion


            #region Body
            
            WANTZ = _lsame.Run(JOBZ, "V");
            LOWER = _lsame.Run(UPLO, "L");
            // *
            INFO = 0;
            if (!(WANTZ || _lsame.Run(JOBZ, "N")))
            {
                INFO =  - 1;
            }
            else
            {
                if (!(LOWER || _lsame.Run(UPLO, "U")))
                {
                    INFO =  - 2;
                }
                else
                {
                    if (N < 0)
                    {
                        INFO =  - 3;
                    }
                    else
                    {
                        if (KD < 0)
                        {
                            INFO =  - 4;
                        }
                        else
                        {
                            if (LDAB < KD + 1)
                            {
                                INFO =  - 6;
                            }
                            else
                            {
                                if (LDZ < 1 || WANTZ && LDZ < N)
                                {
                                    INFO =  - 9;
                                }
                            }
                        }
                    }
                }
            }
            // *
            if (INFO != 0)
            {
                _xerbla.Run("DSBEV ",  - INFO);
                return;
            }
            // *
            // *     Quick return if possible
            // *
            if (N == 0) return;
            // *
            if (N == 1)
            {
                if (LOWER)
                {
                    W[1 + o_w] = AB[1+1 * LDAB + o_ab];
                }
                else
                {
                    W[1 + o_w] = AB[KD + 1+1 * LDAB + o_ab];
                }
                if (WANTZ) Z[1+1 * LDZ + o_z] = ONE;
                return;
            }
            // *
            // *     Get machine constants.
            // *
            SAFMIN = _dlamch.Run("Safe minimum");
            EPS = _dlamch.Run("Precision");
            SMLNUM = SAFMIN / EPS;
            BIGNUM = ONE / SMLNUM;
            RMIN = Math.Sqrt(SMLNUM);
            RMAX = Math.Sqrt(BIGNUM);
            // *
            // *     Scale matrix to allowable range, if necessary.
            // *
            ANRM = _dlansb.Run("M", UPLO, N, KD, AB, offset_ab, LDAB, ref WORK, offset_work);
            ISCALE = 0;
            if (ANRM > ZERO && ANRM < RMIN)
            {
                ISCALE = 1;
                SIGMA = RMIN / ANRM;
            }
            else
            {
                if (ANRM > RMAX)
                {
                    ISCALE = 1;
                    SIGMA = RMAX / ANRM;
                }
            }
            if (ISCALE == 1)
            {
                if (LOWER)
                {
                    _dlascl.Run("B", KD, KD, ONE, SIGMA, N
                                     , N, ref AB, offset_ab, LDAB, ref INFO);
                }
                else
                {
                    _dlascl.Run("Q", KD, KD, ONE, SIGMA, N
                                     , N, ref AB, offset_ab, LDAB, ref INFO);
                }
            }
            // *
            // *     Call DSBTRD to reduce symmetric band matrix to tridiagonal form.
            // *
            INDE = 1;
            INDWRK = INDE + N;
            _dsbtrd.Run(JOBZ, UPLO, N, KD, ref AB, offset_ab, LDAB
                             , ref W, offset_w, ref WORK, INDE + o_work, ref Z, offset_z, LDZ, ref WORK, INDWRK + o_work, ref IINFO);
            // *
            // *     For eigenvalues only, call DSTERF.  For eigenvectors, call SSTEQR.
            // *
            if (!WANTZ)
            {
                _dsterf.Run(N, ref W, offset_w, ref WORK, INDE + o_work, ref INFO);
            }
            else
            {
                _dsteqr.Run(JOBZ, N, ref W, offset_w, ref WORK, INDE + o_work, ref Z, offset_z, LDZ
                                 , ref WORK, INDWRK + o_work, ref INFO);
            }
            // *
            // *     If matrix was scaled, then rescale eigenvalues appropriately.
            // *
            if (ISCALE == 1)
            {
                if (INFO == 0)
                {
                    IMAX = N;
                }
                else
                {
                    IMAX = INFO - 1;
                }
                _dscal.Run(IMAX, ONE / SIGMA, ref W, offset_w, 1);
            }
            // *
            return;
            // *
            // *     End of DSBEV
            // *

            #endregion

        }
    }
}
