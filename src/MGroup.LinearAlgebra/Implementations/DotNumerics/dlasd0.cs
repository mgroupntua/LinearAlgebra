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
    /// Using a divide and conquer approach, DLASD0 computes the singular
    /// value decomposition (SVD) of a real upper bidiagonal N-by-M
    /// matrix B with diagonal D and offdiagonal E, where M = N + SQRE.
    /// The algorithm computes orthogonal matrices U and VT such that
    /// B = U * S * VT. The singular values S are overwritten on D.
    /// 
    /// A related subroutine, DLASDA, computes only the singular values,
    /// and optionally, the singular vectors in compact form.
    /// 
    ///</summary>
    public class DLASD0
    {
    

        #region Dependencies
        
        DLASD1 _dlasd1; DLASDQ _dlasdq; DLASDT _dlasdt; XERBLA _xerbla; 

        #endregion

        public DLASD0(DLASD1 dlasd1, DLASDQ dlasdq, DLASDT dlasdt, XERBLA xerbla)
        {
    

            #region Set Dependencies
            
            _dlasd1 = dlasd1; _dlasdq = dlasdq; _dlasdt = dlasdt; _xerbla = xerbla; 

            #endregion

        }
    
        public DLASD0()
        {
    

            #region Dependencies (Initialization)
            
            var dlamrg = new DLAMRG();
            var lsame = new LSAME();
            var dlamc3 = new DLAMC3();
            var xerbla = new XERBLA();
            var dlapy2 = new DLAPY2();
            var dcopy = new DCOPY();
            var drot = new DROT();
            var dnrm2 = new DNRM2();
            var dlasd5 = new DLASD5();
            var dlas2 = new DLAS2();
            var dlasq5 = new DLASQ5();
            var dlazq4 = new DLAZQ4();
            var ieeeck = new IEEECK();
            var iparmq = new IPARMQ();
            var dscal = new DSCAL();
            var dswap = new DSWAP();
            var dlasdt = new DLASDT();
            var dlamc1 = new DLAMC1(dlamc3);
            var dlamc4 = new DLAMC4(dlamc3);
            var dlamc5 = new DLAMC5(dlamc3);
            var dlamc2 = new DLAMC2(dlamc3, dlamc1, dlamc4, dlamc5);
            var dlamch = new DLAMCH(lsame, dlamc2);
            var dlascl = new DLASCL(lsame, dlamch, xerbla);
            var dlacpy = new DLACPY(lsame);
            var dlaset = new DLASET(lsame);
            var dlasd2 = new DLASD2(dlamch, dlapy2, dcopy, dlacpy, dlamrg, dlaset, drot, xerbla);
            var dgemm = new DGEMM(lsame, xerbla);
            var dlaed6 = new DLAED6(dlamch);
            var dlasd4 = new DLASD4(dlaed6, dlasd5, dlamch);
            var dlasd3 = new DLASD3(dlamc3, dnrm2, dcopy, dgemm, dlacpy, dlascl, dlasd4, xerbla);
            var dlasd1 = new DLASD1(dlamrg, dlascl, dlasd2, dlasd3, xerbla);
            var dlartg = new DLARTG(dlamch);
            var dlasq6 = new DLASQ6(dlamch);
            var dlazq3 = new DLAZQ3(dlasq5, dlasq6, dlazq4, dlamch);
            var dlasrt = new DLASRT(lsame, xerbla);
            var ilaenv = new ILAENV(ieeeck, iparmq);
            var dlasq2 = new DLASQ2(dlazq3, dlasrt, xerbla, dlamch, ilaenv);
            var dlasq1 = new DLASQ1(dcopy, dlas2, dlascl, dlasq2, dlasrt, xerbla, dlamch);
            var dlasr = new DLASR(lsame, xerbla);
            var dlasv2 = new DLASV2(dlamch);
            var dbdsqr = new DBDSQR(lsame, dlamch, dlartg, dlas2, dlasq1, dlasr, dlasv2, drot, dscal, dswap
                                       , xerbla);
            var dlasdq = new DLASDQ(dbdsqr, dlartg, dlasr, dswap, xerbla, lsame);

            #endregion


            #region Set Dependencies
            
            _dlasd1 = dlasd1; _dlasdq = dlasdq; _dlasdt = dlasdt; _xerbla = xerbla; 

            #endregion

        }
        /// <summary>
        /// Purpose
        /// =======
        /// 
        /// Using a divide and conquer approach, DLASD0 computes the singular
        /// value decomposition (SVD) of a real upper bidiagonal N-by-M
        /// matrix B with diagonal D and offdiagonal E, where M = N + SQRE.
        /// The algorithm computes orthogonal matrices U and VT such that
        /// B = U * S * VT. The singular values S are overwritten on D.
        /// 
        /// A related subroutine, DLASDA, computes only the singular values,
        /// and optionally, the singular vectors in compact form.
        /// 
        ///</summary>
        /// <param name="N">
        /// (input) INTEGER
        /// On entry, the row dimension of the upper bidiagonal matrix.
        /// This is also the dimension of the main diagonal array D.
        ///</param>
        /// <param name="SQRE">
        /// (input) INTEGER
        /// Specifies the column dimension of the bidiagonal matrix.
        /// = 0: The bidiagonal matrix has column dimension M = N;
        /// = 1: The bidiagonal matrix has column dimension M = N+1;
        ///</param>
        /// <param name="D">
        /// (input/output) DOUBLE PRECISION array, dimension (N)
        /// On entry D contains the main diagonal of the bidiagonal
        /// matrix.
        /// On exit D, if INFO = 0, contains its singular values.
        ///</param>
        /// <param name="E">
        /// (input) DOUBLE PRECISION array, dimension (M-1)
        /// Contains the subdiagonal entries of the bidiagonal matrix.
        /// On exit, E has been destroyed.
        ///</param>
        /// <param name="U">
        /// (output) DOUBLE PRECISION array, dimension at least (LDQ, N)
        /// On exit, U contains the left singular vectors.
        ///</param>
        /// <param name="LDU">
        /// (input) INTEGER
        /// On entry, leading dimension of U.
        ///</param>
        /// <param name="VT">
        /// (output) DOUBLE PRECISION array, dimension at least (LDVT, M)
        /// On exit, VT' contains the right singular vectors.
        ///</param>
        /// <param name="LDVT">
        /// (input) INTEGER
        /// On entry, leading dimension of VT.
        ///</param>
        /// <param name="SMLSIZ">
        /// (input) INTEGER
        /// On entry, maximum size of the subproblems at the
        /// bottom of the computation tree.
        ///</param>
        /// <param name="IWORK">
        /// (workspace) INTEGER work array.
        /// Dimension must be at least (8 * N)
        ///</param>
        /// <param name="WORK">
        /// (workspace) DOUBLE PRECISION work array.
        /// Dimension must be at least (3 * M**2 + 2 * M)
        ///</param>
        /// <param name="INFO">
        /// (output) INTEGER
        /// = 0:  successful exit.
        /// .LT. 0:  if INFO = -i, the i-th argument had an illegal value.
        /// .GT. 0:  if INFO = 1, an singular value did not converge
        ///</param>
        public void Run(int N, int SQRE, ref double[] D, int offset_d, ref double[] E, int offset_e, ref double[] U, int offset_u, int LDU
                         , ref double[] VT, int offset_vt, int LDVT, int SMLSIZ, ref int[] IWORK, int offset_iwork, ref double[] WORK, int offset_work, ref int INFO)
        {

            #region Variables
            
            var I = 0; var I1 = 0; var IC = 0; var IDXQ = 0; var IDXQC = 0; var IM1 = 0; var INODE = 0; var ITEMP = 0; 
            var IWK = 0;var J = 0; var LF = 0; var LL = 0; var LVL = 0; var M = 0; var NCC = 0; var ND = 0; var NDB1 = 0; 
            var NDIML = 0;var NDIMR = 0; var NL = 0; var NLF = 0; var NLP1 = 0; var NLVL = 0; var NR = 0; var NRF = 0; 
            var NRP1 = 0;var SQREI = 0; double ALPHA = 0; double BETA = 0; 

            #endregion


            #region Array Index Correction
            
             var o_d = -1 + offset_d;  var o_e = -1 + offset_e;  var o_u = -1 - LDU + offset_u;  var o_vt = -1 - LDVT + offset_vt; 
             var o_iwork = -1 + offset_iwork; var o_work = -1 + offset_work; 

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
            // *  Using a divide and conquer approach, DLASD0 computes the singular
            // *  value decomposition (SVD) of a real upper bidiagonal N-by-M
            // *  matrix B with diagonal D and offdiagonal E, where M = N + SQRE.
            // *  The algorithm computes orthogonal matrices U and VT such that
            // *  B = U * S * VT. The singular values S are overwritten on D.
            // *
            // *  A related subroutine, DLASDA, computes only the singular values,
            // *  and optionally, the singular vectors in compact form.
            // *
            // *  Arguments
            // *  =========
            // *
            // *  N      (input) INTEGER
            // *         On entry, the row dimension of the upper bidiagonal matrix.
            // *         This is also the dimension of the main diagonal array D.
            // *
            // *  SQRE   (input) INTEGER
            // *         Specifies the column dimension of the bidiagonal matrix.
            // *         = 0: The bidiagonal matrix has column dimension M = N;
            // *         = 1: The bidiagonal matrix has column dimension M = N+1;
            // *
            // *  D      (input/output) DOUBLE PRECISION array, dimension (N)
            // *         On entry D contains the main diagonal of the bidiagonal
            // *         matrix.
            // *         On exit D, if INFO = 0, contains its singular values.
            // *
            // *  E      (input) DOUBLE PRECISION array, dimension (M-1)
            // *         Contains the subdiagonal entries of the bidiagonal matrix.
            // *         On exit, E has been destroyed.
            // *
            // *  U      (output) DOUBLE PRECISION array, dimension at least (LDQ, N)
            // *         On exit, U contains the left singular vectors.
            // *
            // *  LDU    (input) INTEGER
            // *         On entry, leading dimension of U.
            // *
            // *  VT     (output) DOUBLE PRECISION array, dimension at least (LDVT, M)
            // *         On exit, VT' contains the right singular vectors.
            // *
            // *  LDVT   (input) INTEGER
            // *         On entry, leading dimension of VT.
            // *
            // *  SMLSIZ (input) INTEGER
            // *         On entry, maximum size of the subproblems at the
            // *         bottom of the computation tree.
            // *
            // *  IWORK  (workspace) INTEGER work array.
            // *         Dimension must be at least (8 * N)
            // *
            // *  WORK   (workspace) DOUBLE PRECISION work array.
            // *         Dimension must be at least (3 * M**2 + 2 * M)
            // *
            // *  INFO   (output) INTEGER
            // *          = 0:  successful exit.
            // *          < 0:  if INFO = -i, the i-th argument had an illegal value.
            // *          > 0:  if INFO = 1, an singular value did not converge
            // *
            // *  Further Details
            // *  ===============
            // *
            // *  Based on contributions by
            // *     Ming Gu and Huan Ren, Computer Science Division, University of
            // *     California at Berkeley, USA
            // *
            // *  =====================================================================
            // *
            // *     .. Local Scalars ..
            // *     ..
            // *     .. External Subroutines ..
            // *     ..
            // *     .. Executable Statements ..
            // *
            // *     Test the input parameters.
            // *

            #endregion


            #region Body
            
            INFO = 0;
            // *
            if (N < 0)
            {
                INFO =  - 1;
            }
            else
            {
                if (SQRE < 0 || SQRE > 1)
                {
                    INFO =  - 2;
                }
            }
            // *
            M = N + SQRE;
            // *
            if (LDU < N)
            {
                INFO =  - 6;
            }
            else
            {
                if (LDVT < M)
                {
                    INFO =  - 8;
                }
                else
                {
                    if (SMLSIZ < 3)
                    {
                        INFO =  - 9;
                    }
                }
            }
            if (INFO != 0)
            {
                _xerbla.Run("DLASD0",  - INFO);
                return;
            }
            // *
            // *     If the input matrix is too small, call DLASDQ to find the SVD.
            // *
            if (N <= SMLSIZ)
            {
                _dlasdq.Run("U", SQRE, N, M, N, 0
                                 , ref D, offset_d, ref E, offset_e, ref VT, offset_vt, LDVT, ref U, offset_u, LDU
                                 , ref U, offset_u, LDU, ref WORK, offset_work, ref INFO);
                return;
            }
            // *
            // *     Set up the computation tree.
            // *
            INODE = 1;
            NDIML = INODE + N;
            NDIMR = NDIML + N;
            IDXQ = NDIMR + N;
            IWK = IDXQ + N;
            _dlasdt.Run(N, ref NLVL, ref ND, ref IWORK, INODE + o_iwork, ref IWORK, NDIML + o_iwork, ref IWORK, NDIMR + o_iwork
                             , SMLSIZ);
            // *
            // *     For the nodes on bottom level of the tree, solve
            // *     their subproblems by DLASDQ.
            // *
            NDB1 = (ND + 1) / 2;
            NCC = 0;
            for (I = NDB1; I <= ND; I++)
            {
                // *
                // *     IC : center row of each node
                // *     NL : number of rows of left  subproblem
                // *     NR : number of rows of right subproblem
                // *     NLF: starting row of the left   subproblem
                // *     NRF: starting row of the right  subproblem
                // *
                I1 = I - 1;
                IC = IWORK[INODE + I1 + o_iwork];
                NL = IWORK[NDIML + I1 + o_iwork];
                NLP1 = NL + 1;
                NR = IWORK[NDIMR + I1 + o_iwork];
                NRP1 = NR + 1;
                NLF = IC - NL;
                NRF = IC + 1;
                SQREI = 1;
                _dlasdq.Run("U", SQREI, NL, NLP1, NL, NCC
                                 , ref D, NLF + o_d, ref E, NLF + o_e, ref VT, NLF+NLF * LDVT + o_vt, LDVT, ref U, NLF+NLF * LDU + o_u, LDU
                                 , ref U, NLF+NLF * LDU + o_u, LDU, ref WORK, offset_work, ref INFO);
                if (INFO != 0)
                {
                    return;
                }
                ITEMP = IDXQ + NLF - 2;
                for (J = 1; J <= NL; J++)
                {
                    IWORK[ITEMP + J + o_iwork] = J;
                }
                if (I == ND)
                {
                    SQREI = SQRE;
                }
                else
                {
                    SQREI = 1;
                }
                NRP1 = NR + SQREI;
                _dlasdq.Run("U", SQREI, NR, NRP1, NR, NCC
                                 , ref D, NRF + o_d, ref E, NRF + o_e, ref VT, NRF+NRF * LDVT + o_vt, LDVT, ref U, NRF+NRF * LDU + o_u, LDU
                                 , ref U, NRF+NRF * LDU + o_u, LDU, ref WORK, offset_work, ref INFO);
                if (INFO != 0)
                {
                    return;
                }
                ITEMP = IDXQ + IC;
                for (J = 1; J <= NR; J++)
                {
                    IWORK[ITEMP + J - 1 + o_iwork] = J;
                }
            }
            // *
            // *     Now conquer each subproblem bottom-up.
            // *
            for (LVL = NLVL; LVL >= 1; LVL +=  - 1)
            {
                // *
                // *        Find the first node LF and last node LL on the
                // *        current level LVL.
                // *
                if (LVL == 1)
                {
                    LF = 1;
                    LL = 1;
                }
                else
                {
                    LF = (int)Math.Pow(2, LVL - 1);
                    LL = 2 * LF - 1;
                }
                for (I = LF; I <= LL; I++)
                {
                    IM1 = I - 1;
                    IC = IWORK[INODE + IM1 + o_iwork];
                    NL = IWORK[NDIML + IM1 + o_iwork];
                    NR = IWORK[NDIMR + IM1 + o_iwork];
                    NLF = IC - NL;
                    if (SQRE == 0 && I == LL)
                    {
                        SQREI = SQRE;
                    }
                    else
                    {
                        SQREI = 1;
                    }
                    IDXQC = IDXQ + NLF - 1;
                    ALPHA = D[IC + o_d];
                    BETA = E[IC + o_e];
                    _dlasd1.Run(NL, NR, SQREI, ref D, NLF + o_d, ref ALPHA, ref BETA
                                     , ref U, NLF+NLF * LDU + o_u, LDU, ref VT, NLF+NLF * LDVT + o_vt, LDVT, ref IWORK, IDXQC + o_iwork, ref IWORK, IWK + o_iwork
                                     , ref WORK, offset_work, ref INFO);
                    if (INFO != 0)
                    {
                        return;
                    }
                }
            }
            // *
            return;
            // *
            // *     End of DLASD0
            // *

            #endregion

        }
    }
}
