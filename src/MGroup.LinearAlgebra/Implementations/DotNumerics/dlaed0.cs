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
    /// DLAED0 computes all eigenvalues and corresponding eigenvectors of a
    /// symmetric tridiagonal matrix using the divide and conquer method.
    /// 
    ///</summary>
    public class DLAED0
    {
    

        #region Dependencies
        
        DCOPY _dcopy; DGEMM _dgemm; DLACPY _dlacpy; DLAED1 _dlaed1; DLAED7 _dlaed7; DSTEQR _dsteqr; XERBLA _xerbla; 
        ILAENV _ilaenv;

        #endregion


        #region Variables
        
        const double ZERO = 0.0E0; const double ONE = 1.0E0; const double TWO = 2.0E0; 

        #endregion

        public DLAED0(DCOPY dcopy, DGEMM dgemm, DLACPY dlacpy, DLAED1 dlaed1, DLAED7 dlaed7, DSTEQR dsteqr, XERBLA xerbla, ILAENV ilaenv)
        {
    

            #region Set Dependencies
            
            _dcopy = dcopy; _dgemm = dgemm; _dlacpy = dlacpy; _dlaed1 = dlaed1; _dlaed7 = dlaed7; 
            _dsteqr = dsteqr;_xerbla = xerbla; _ilaenv = ilaenv; 

            #endregion

        }
    
        public DLAED0()
        {
    

            #region Dependencies (Initialization)
            
            var dcopy = new DCOPY();
            var lsame = new LSAME();
            var xerbla = new XERBLA();
            var idamax = new IDAMAX();
            var dlamc3 = new DLAMC3();
            var dlapy2 = new DLAPY2();
            var dlamrg = new DLAMRG();
            var drot = new DROT();
            var dscal = new DSCAL();
            var dnrm2 = new DNRM2();
            var dlaed5 = new DLAED5();
            var dlassq = new DLASSQ();
            var dlae2 = new DLAE2();
            var dlaev2 = new DLAEV2();
            var dswap = new DSWAP();
            var ieeeck = new IEEECK();
            var iparmq = new IPARMQ();
            var dgemm = new DGEMM(lsame, xerbla);
            var dlacpy = new DLACPY(lsame);
            var dlamc1 = new DLAMC1(dlamc3);
            var dlamc4 = new DLAMC4(dlamc3);
            var dlamc5 = new DLAMC5(dlamc3);
            var dlamc2 = new DLAMC2(dlamc3, dlamc1, dlamc4, dlamc5);
            var dlamch = new DLAMCH(lsame, dlamc2);
            var dlaed2 = new DLAED2(idamax, dlamch, dlapy2, dcopy, dlacpy, dlamrg, drot, dscal, xerbla);
            var dlaed6 = new DLAED6(dlamch);
            var dlaed4 = new DLAED4(dlamch, dlaed5, dlaed6);
            var dlaset = new DLASET(lsame);
            var dlaed3 = new DLAED3(dlamc3, dnrm2, dcopy, dgemm, dlacpy, dlaed4, dlaset, xerbla);
            var dlaed1 = new DLAED1(dcopy, dlaed2, dlaed3, dlamrg, xerbla);
            var dlaed8 = new DLAED8(idamax, dlamch, dlapy2, dcopy, dlacpy, dlamrg, drot, dscal, xerbla);
            var dlaed9 = new DLAED9(dlamc3, dnrm2, dcopy, dlaed4, xerbla);
            var dgemv = new DGEMV(lsame, xerbla);
            var dlaeda = new DLAEDA(dcopy, dgemv, drot, xerbla);
            var dlaed7 = new DLAED7(dgemm, dlaed8, dlaed9, dlaeda, dlamrg, xerbla);
            var dlanst = new DLANST(lsame, dlassq);
            var dlartg = new DLARTG(dlamch);
            var dlascl = new DLASCL(lsame, dlamch, xerbla);
            var dlasr = new DLASR(lsame, xerbla);
            var dlasrt = new DLASRT(lsame, xerbla);
            var dsteqr = new DSTEQR(lsame, dlamch, dlanst, dlapy2, dlae2, dlaev2, dlartg, dlascl, dlaset, dlasr
                                       , dlasrt, dswap, xerbla);
            var ilaenv = new ILAENV(ieeeck, iparmq);

            #endregion


            #region Set Dependencies
            
            _dcopy = dcopy; _dgemm = dgemm; _dlacpy = dlacpy; _dlaed1 = dlaed1; _dlaed7 = dlaed7; 
            _dsteqr = dsteqr;_xerbla = xerbla; _ilaenv = ilaenv; 

            #endregion

        }
        /// <summary>
        /// Purpose
        /// =======
        /// 
        /// DLAED0 computes all eigenvalues and corresponding eigenvectors of a
        /// symmetric tridiagonal matrix using the divide and conquer method.
        /// 
        ///</summary>
        /// <param name="ICOMPQ">
        /// (input) INTEGER
        /// = 0:  Compute eigenvalues only.
        /// = 1:  Compute eigenvectors of original dense symmetric matrix
        /// also.  On entry, Q contains the orthogonal matrix used
        /// to reduce the original matrix to tridiagonal form.
        /// = 2:  Compute eigenvalues and eigenvectors of tridiagonal
        /// matrix.
        ///</param>
        /// <param name="QSIZ">
        /// (input) INTEGER
        /// The dimension of the orthogonal matrix used to reduce
        /// the full matrix to tridiagonal form.  QSIZ .GE. N if ICOMPQ = 1.
        ///</param>
        /// <param name="N">
        /// (input) INTEGER
        /// The dimension of the symmetric tridiagonal matrix.  N .GE. 0.
        ///</param>
        /// <param name="D">
        /// (input/output) DOUBLE PRECISION array, dimension (N)
        /// On entry, the main diagonal of the tridiagonal matrix.
        /// On exit, its eigenvalues.
        ///</param>
        /// <param name="E">
        /// (input) DOUBLE PRECISION array, dimension (N-1)
        /// The off-diagonal elements of the tridiagonal matrix.
        /// On exit, E has been destroyed.
        ///</param>
        /// <param name="Q">
        /// (input/output) DOUBLE PRECISION array, dimension (LDQ, N)
        /// On entry, Q must contain an N-by-N orthogonal matrix.
        /// If ICOMPQ = 0    Q is not referenced.
        /// If ICOMPQ = 1    On entry, Q is a subset of the columns of the
        /// orthogonal matrix used to reduce the full
        /// matrix to tridiagonal form corresponding to
        /// the subset of the full matrix which is being
        /// decomposed at this time.
        /// If ICOMPQ = 2    On entry, Q will be the identity matrix.
        /// On exit, Q contains the eigenvectors of the
        /// tridiagonal matrix.
        ///</param>
        /// <param name="LDQ">
        /// (input) INTEGER
        /// The leading dimension of the array Q.  If eigenvectors are
        /// desired, then  LDQ .GE. max(1,N).  In any case,  LDQ .GE. 1.
        ///</param>
        /// <param name="QSTORE">
        /// (workspace) DOUBLE PRECISION array, dimension (LDQS, N)
        /// Referenced only when ICOMPQ = 1.  Used to store parts of
        /// the eigenvector matrix when the updating matrix multiplies
        /// take place.
        ///</param>
        /// <param name="LDQS">
        /// (input) INTEGER
        /// The leading dimension of the array QSTORE.  If ICOMPQ = 1,
        /// then  LDQS .GE. max(1,N).  In any case,  LDQS .GE. 1.
        ///</param>
        /// <param name="WORK">
        /// (workspace) DOUBLE PRECISION array,
        /// If ICOMPQ = 0 or 1, the dimension of WORK must be at least
        /// 1 + 3*N + 2*N*lg N + 2*N**2
        /// ( lg( N ) = smallest integer k
        /// such that 2^k .GE. N )
        /// If ICOMPQ = 2, the dimension of WORK must be at least
        /// 4*N + N**2.
        ///</param>
        /// <param name="IWORK">
        /// (workspace) INTEGER array,
        /// If ICOMPQ = 0 or 1, the dimension of IWORK must be at least
        /// 6 + 6*N + 5*N*lg N.
        /// ( lg( N ) = smallest integer k
        /// such that 2^k .GE. N )
        /// If ICOMPQ = 2, the dimension of IWORK must be at least
        /// 3 + 5*N.
        ///</param>
        /// <param name="INFO">
        /// (output) INTEGER
        /// = 0:  successful exit.
        /// .LT. 0:  if INFO = -i, the i-th argument had an illegal value.
        /// .GT. 0:  The algorithm failed to compute an eigenvalue while
        /// working on the submatrix lying in rows and columns
        /// INFO/(N+1) through mod(INFO,N+1).
        ///</param>
        public void Run(int ICOMPQ, int QSIZ, int N, ref double[] D, int offset_d, ref double[] E, int offset_e, ref double[] Q, int offset_q
                         , int LDQ, ref double[] QSTORE, int offset_qstore, int LDQS, ref double[] WORK, int offset_work, ref int[] IWORK, int offset_iwork, ref int INFO)
        {

            #region Variables
            
            var CURLVL = 0; var CURPRB = 0; var CURR = 0; var I = 0; var IGIVCL = 0; var IGIVNM = 0; var IGIVPT = 0; 
            var INDXQ = 0;var IPERM = 0; var IPRMPT = 0; var IQ = 0; var IQPTR = 0; var IWREM = 0; var J = 0; var K = 0; 
            var LGN = 0;var MATSIZ = 0; var MSD2 = 0; var SMLSIZ = 0; var SMM1 = 0; var SPM1 = 0; var SPM2 = 0; var SUBMAT = 0; 
            var SUBPBS = 0;var TLVLS = 0; double TEMP = 0; 

            #endregion


            #region Array Index Correction
            
             var o_d = -1 + offset_d;  var o_e = -1 + offset_e;  var o_q = -1 - LDQ + offset_q; 
             var o_qstore = -1 - LDQS + offset_qstore; var o_work = -1 + offset_work;  var o_iwork = -1 + offset_iwork; 

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
            // *  DLAED0 computes all eigenvalues and corresponding eigenvectors of a
            // *  symmetric tridiagonal matrix using the divide and conquer method.
            // *
            // *  Arguments
            // *  =========
            // *
            // *  ICOMPQ  (input) INTEGER
            // *          = 0:  Compute eigenvalues only.
            // *          = 1:  Compute eigenvectors of original dense symmetric matrix
            // *                also.  On entry, Q contains the orthogonal matrix used
            // *                to reduce the original matrix to tridiagonal form.
            // *          = 2:  Compute eigenvalues and eigenvectors of tridiagonal
            // *                matrix.
            // *
            // *  QSIZ   (input) INTEGER
            // *         The dimension of the orthogonal matrix used to reduce
            // *         the full matrix to tridiagonal form.  QSIZ >= N if ICOMPQ = 1.
            // *
            // *  N      (input) INTEGER
            // *         The dimension of the symmetric tridiagonal matrix.  N >= 0.
            // *
            // *  D      (input/output) DOUBLE PRECISION array, dimension (N)
            // *         On entry, the main diagonal of the tridiagonal matrix.
            // *         On exit, its eigenvalues.
            // *
            // *  E      (input) DOUBLE PRECISION array, dimension (N-1)
            // *         The off-diagonal elements of the tridiagonal matrix.
            // *         On exit, E has been destroyed.
            // *
            // *  Q      (input/output) DOUBLE PRECISION array, dimension (LDQ, N)
            // *         On entry, Q must contain an N-by-N orthogonal matrix.
            // *         If ICOMPQ = 0    Q is not referenced.
            // *         If ICOMPQ = 1    On entry, Q is a subset of the columns of the
            // *                          orthogonal matrix used to reduce the full
            // *                          matrix to tridiagonal form corresponding to
            // *                          the subset of the full matrix which is being
            // *                          decomposed at this time.
            // *         If ICOMPQ = 2    On entry, Q will be the identity matrix.
            // *                          On exit, Q contains the eigenvectors of the
            // *                          tridiagonal matrix.
            // *
            // *  LDQ    (input) INTEGER
            // *         The leading dimension of the array Q.  If eigenvectors are
            // *         desired, then  LDQ >= max(1,N).  In any case,  LDQ >= 1.
            // *
            // *  QSTORE (workspace) DOUBLE PRECISION array, dimension (LDQS, N)
            // *         Referenced only when ICOMPQ = 1.  Used to store parts of
            // *         the eigenvector matrix when the updating matrix multiplies
            // *         take place.
            // *
            // *  LDQS   (input) INTEGER
            // *         The leading dimension of the array QSTORE.  If ICOMPQ = 1,
            // *         then  LDQS >= max(1,N).  In any case,  LDQS >= 1.
            // *
            // *  WORK   (workspace) DOUBLE PRECISION array,
            // *         If ICOMPQ = 0 or 1, the dimension of WORK must be at least
            // *                     1 + 3*N + 2*N*lg N + 2*N**2
            // *                     ( lg( N ) = smallest integer k
            // *                                 such that 2^k >= N )
            // *         If ICOMPQ = 2, the dimension of WORK must be at least
            // *                     4*N + N**2.
            // *
            // *  IWORK  (workspace) INTEGER array,
            // *         If ICOMPQ = 0 or 1, the dimension of IWORK must be at least
            // *                        6 + 6*N + 5*N*lg N.
            // *                        ( lg( N ) = smallest integer k
            // *                                    such that 2^k >= N )
            // *         If ICOMPQ = 2, the dimension of IWORK must be at least
            // *                        3 + 5*N.
            // *
            // *  INFO   (output) INTEGER
            // *          = 0:  successful exit.
            // *          < 0:  if INFO = -i, the i-th argument had an illegal value.
            // *          > 0:  The algorithm failed to compute an eigenvalue while
            // *                working on the submatrix lying in rows and columns
            // *                INFO/(N+1) through mod(INFO,N+1).
            // *
            // *  Further Details
            // *  ===============
            // *
            // *  Based on contributions by
            // *     Jeff Rutter, Computer Science Division, University of California
            // *     at Berkeley, USA
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
            //      INTRINSIC          ABS, DBLE, INT, LOG, MAX;
            // *     ..
            // *     .. Executable Statements ..
            // *
            // *     Test the input parameters.
            // *

            #endregion


            #region Body
            
            INFO = 0;
            // *
            if (ICOMPQ < 0 || ICOMPQ > 2)
            {
                INFO =  - 1;
            }
            else
            {
                if (ICOMPQ == 1 && QSIZ < Math.Max(0, N))
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
                        if (LDQ < Math.Max(1, N))
                        {
                            INFO =  - 7;
                        }
                        else
                        {
                            if (LDQS < Math.Max(1, N))
                            {
                                INFO =  - 9;
                            }
                        }
                    }
                }
            }
            if (INFO != 0)
            {
                _xerbla.Run("DLAED0",  - INFO);
                return;
            }
            // *
            // *     Quick return if possible
            // *
            if (N == 0) return;
            // *
            SMLSIZ = _ilaenv.Run(9, "DLAED0", " ", 0, 0, 0, 0);
            // *
            // *     Determine the size and placement of the submatrices, and save in
            // *     the leading elements of IWORK.
            // *
            IWORK[1 + o_iwork] = N;
            SUBPBS = 1;
            TLVLS = 0;
        LABEL10:;
            if (IWORK[SUBPBS + o_iwork] > SMLSIZ)
            {
                for (J = SUBPBS; J >= 1; J +=  - 1)
                {
                    IWORK[2 * J + o_iwork] = (IWORK[J + o_iwork] + 1) / 2;
                    IWORK[2 * J - 1 + o_iwork] = IWORK[J + o_iwork] / 2;
                }
                TLVLS += 1;
                SUBPBS *= 2;
                goto LABEL10;
            }
            for (J = 2; J <= SUBPBS; J++)
            {
                IWORK[J + o_iwork] += IWORK[J - 1 + o_iwork];
            }
            // *
            // *     Divide the matrix into SUBPBS submatrices of size at most SMLSIZ+1
            // *     using rank-1 modifications (cuts).
            // *
            SPM1 = SUBPBS - 1;
            for (I = 1; I <= SPM1; I++)
            {
                SUBMAT = IWORK[I + o_iwork] + 1;
                SMM1 = SUBMAT - 1;
                D[SMM1 + o_d] -= Math.Abs(E[SMM1 + o_e]);
                D[SUBMAT + o_d] -= Math.Abs(E[SMM1 + o_e]);
            }
            // *
            INDXQ = 4 * N + 3;
            if (ICOMPQ != 2)
            {
                // *
                // *        Set up workspaces for eigenvalues only/accumulate new vectors
                // *        routine
                // *
                TEMP = Math.Log(Convert.ToDouble(N)) / Math.Log(TWO);
                LGN = Convert.ToInt32(Math.Truncate(TEMP));
                if (Math.Pow(2,LGN) < N) LGN += 1;
                if (Math.Pow(2,LGN) < N) LGN += 1;
                IPRMPT = INDXQ + N + 1;
                IPERM = IPRMPT + N * LGN;
                IQPTR = IPERM + N * LGN;
                IGIVPT = IQPTR + N + 2;
                IGIVCL = IGIVPT + N * LGN;
                // *
                IGIVNM = 1;
                IQ = IGIVNM + 2 * N * LGN;
                IWREM = IQ + (int)Math.Pow(N, 2) + 1;
                // *
                // *        Initialize pointers
                // *
                for (I = 0; I <= SUBPBS; I++)
                {
                    IWORK[IPRMPT + I + o_iwork] = 1;
                    IWORK[IGIVPT + I + o_iwork] = 1;
                }
                IWORK[IQPTR + o_iwork] = 1;
            }
            // *
            // *     Solve each submatrix eigenproblem at the bottom of the divide and
            // *     conquer tree.
            // *
            CURR = 0;
            for (I = 0; I <= SPM1; I++)
            {
                if (I == 0)
                {
                    SUBMAT = 1;
                    MATSIZ = IWORK[1 + o_iwork];
                }
                else
                {
                    SUBMAT = IWORK[I + o_iwork] + 1;
                    MATSIZ = IWORK[I + 1 + o_iwork] - IWORK[I + o_iwork];
                }
                if (ICOMPQ == 2)
                {
                    _dsteqr.Run("I", MATSIZ, ref D, SUBMAT + o_d, ref E, SUBMAT + o_e, ref Q, SUBMAT+SUBMAT * LDQ + o_q, LDQ
                                     , ref WORK, offset_work, ref INFO);
                    if (INFO != 0) goto LABEL130;
                }
                else
                {
                    _dsteqr.Run("I", MATSIZ, ref D, SUBMAT + o_d, ref E, SUBMAT + o_e, ref WORK, IQ - 1 + IWORK[IQPTR + CURR + o_iwork] + o_work, MATSIZ
                                     , ref WORK, offset_work, ref INFO);
                    if (INFO != 0) goto LABEL130;
                    if (ICOMPQ == 1)
                    {
                        _dgemm.Run("N", "N", QSIZ, MATSIZ, MATSIZ, ONE
                                        , Q, 1+SUBMAT * LDQ + o_q, LDQ, WORK, IQ - 1 + IWORK[IQPTR + CURR + o_iwork] + o_work, MATSIZ, ZERO, ref QSTORE, 1+SUBMAT * LDQS + o_qstore
                                        , LDQS);
                    }
                    IWORK[IQPTR + CURR + 1 + o_iwork] = IWORK[IQPTR + CURR + o_iwork] + (int)Math.Pow(MATSIZ, 2);
                    CURR += 1;
                }
                K = 1;
                for (J = SUBMAT; J <= IWORK[I + 1 + o_iwork]; J++)
                {
                    IWORK[INDXQ + J + o_iwork] = K;
                    K += 1;
                }
            }
            // *
            // *     Successively merge eigensystems of adjacent submatrices
            // *     into eigensystem for the corresponding larger matrix.
            // *
            // *     while ( SUBPBS > 1 )
            // *
            CURLVL = 1;
        LABEL80:;
            if (SUBPBS > 1)
            {
                SPM2 = SUBPBS - 2;
                for (I = 0; I <= SPM2; I += 2)
                {
                    if (I == 0)
                    {
                        SUBMAT = 1;
                        MATSIZ = IWORK[2 + o_iwork];
                        MSD2 = IWORK[1 + o_iwork];
                        CURPRB = 0;
                    }
                    else
                    {
                        SUBMAT = IWORK[I + o_iwork] + 1;
                        MATSIZ = IWORK[I + 2 + o_iwork] - IWORK[I + o_iwork];
                        MSD2 = MATSIZ / 2;
                        CURPRB += 1;
                    }
                    // *
                    // *     Merge lower order eigensystems (of size MSD2 and MATSIZ - MSD2)
                    // *     into an eigensystem of size MATSIZ.
                    // *     DLAED1 is used only for the full eigensystem of a tridiagonal
                    // *     matrix.
                    // *     DLAED7 handles the cases in which eigenvalues only or eigenvalues
                    // *     and eigenvectors of a full symmetric matrix (which was reduced to
                    // *     tridiagonal form) are desired.
                    // *
                    if (ICOMPQ == 2)
                    {
                        _dlaed1.Run(MATSIZ, ref D, SUBMAT + o_d, ref Q, SUBMAT+SUBMAT * LDQ + o_q, LDQ, ref IWORK, INDXQ + SUBMAT + o_iwork, ref E[SUBMAT + MSD2 - 1 + o_e]
                                         , MSD2, ref WORK, offset_work, ref IWORK, SUBPBS + 1 + o_iwork, ref INFO);
                    }
                    else
                    {
                        _dlaed7.Run(ICOMPQ, MATSIZ, QSIZ, TLVLS, CURLVL, CURPRB
                                         , ref D, SUBMAT + o_d, ref QSTORE, 1+SUBMAT * LDQS + o_qstore, LDQS, ref IWORK, INDXQ + SUBMAT + o_iwork, ref E[SUBMAT + MSD2 - 1 + o_e], MSD2
                                         , ref WORK, IQ + o_work, ref IWORK, IQPTR + o_iwork, ref IWORK, IPRMPT + o_iwork, ref IWORK, IPERM + o_iwork, ref IWORK, IGIVPT + o_iwork, ref IWORK, IGIVCL + o_iwork
                                         , ref WORK, IGIVNM + o_work, ref WORK, IWREM + o_work, ref IWORK, SUBPBS + 1 + o_iwork, ref INFO);
                    }
                    if (INFO != 0) goto LABEL130;
                    IWORK[I / 2 + 1 + o_iwork] = IWORK[I + 2 + o_iwork];
                }
                SUBPBS /= 2;
                CURLVL += 1;
                goto LABEL80;
            }
            // *
            // *     end while
            // *
            // *     Re-merge the eigenvalues/vectors which were deflated at the final
            // *     merge step.
            // *
            if (ICOMPQ == 1)
            {
                for (I = 1; I <= N; I++)
                {
                    J = IWORK[INDXQ + I + o_iwork];
                    WORK[I + o_work] = D[J + o_d];
                    _dcopy.Run(QSIZ, QSTORE, 1+J * LDQS + o_qstore, 1, ref Q, 1+I * LDQ + o_q, 1);
                }
                _dcopy.Run(N, WORK, offset_work, 1, ref D, offset_d, 1);
            }
            else
            {
                if (ICOMPQ == 2)
                {
                    for (I = 1; I <= N; I++)
                    {
                        J = IWORK[INDXQ + I + o_iwork];
                        WORK[I + o_work] = D[J + o_d];
                        _dcopy.Run(N, Q, 1+J * LDQ + o_q, 1, ref WORK, N * I + 1 + o_work, 1);
                    }
                    _dcopy.Run(N, WORK, offset_work, 1, ref D, offset_d, 1);
                    _dlacpy.Run("A", N, N, WORK, N + 1 + o_work, N, ref Q, offset_q
                                     , LDQ);
                }
                else
                {
                    for (I = 1; I <= N; I++)
                    {
                        J = IWORK[INDXQ + I + o_iwork];
                        WORK[I + o_work] = D[J + o_d];
                    }
                    _dcopy.Run(N, WORK, offset_work, 1, ref D, offset_d, 1);
                }
            }
            goto LABEL140;
            // *
        LABEL130:;
            INFO = SUBMAT * (N + 1) + SUBMAT + MATSIZ - 1;
            // *
        LABEL140:;
            return;
            // *
            // *     End of DLAED0
            // *

            #endregion

        }
    }
}
