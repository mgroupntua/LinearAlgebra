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
    /// DGESDD computes the singular value decomposition (SVD) of a real
    /// M-by-N matrix A, optionally computing the left and right singular
    /// vectors.  If singular vectors are desired, it uses a
    /// divide-and-conquer algorithm.
    /// 
    /// The SVD is written
    /// 
    /// A = U * SIGMA * transpose(V)
    /// 
    /// where SIGMA is an M-by-N matrix which is zero except for its
    /// min(m,n) diagonal elements, U is an M-by-M orthogonal matrix, and
    /// V is an N-by-N orthogonal matrix.  The diagonal elements of SIGMA
    /// are the singular values of A; they are real and non-negative, and
    /// are returned in descending order.  The first min(m,n) columns of
    /// U and V are the left and right singular vectors of A.
    /// 
    /// Note that the routine returns VT = V**T, not V.
    /// 
    /// The divide and conquer algorithm makes very mild assumptions about
    /// floating point arithmetic. It will work on machines with a guard
    /// digit in add/subtract, or on those binary machines without guard
    /// digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
    /// Cray-2. It could conceivably fail on hexadecimal or decimal machines
    /// without guard digits, but we know of none.
    /// 
    ///</summary>
    public class DGESDD
    {
    

        #region Dependencies
        
        DBDSDC _dbdsdc; DGEBRD _dgebrd; DGELQF _dgelqf; DGEMM _dgemm; DGEQRF _dgeqrf; DLACPY _dlacpy; DLASCL _dlascl; 
        DLASET _dlaset;DORGBR _dorgbr; DORGLQ _dorglq; DORGQR _dorgqr; DORMBR _dormbr; XERBLA _xerbla; DLAMCH _dlamch; 
        DLANGE _dlange;ILAENV _ilaenv; LSAME _lsame; 

        #endregion


        #region Variables
        
        const double ZERO = 0.0E0; const double ONE = 1.0E0; int[] IDUM = new int[1]; double[] DUM = new double[1]; 

        #endregion

        public DGESDD(DBDSDC dbdsdc, DGEBRD dgebrd, DGELQF dgelqf, DGEMM dgemm, DGEQRF dgeqrf, DLACPY dlacpy, DLASCL dlascl, DLASET dlaset, DORGBR dorgbr, DORGLQ dorglq
                      , DORGQR dorgqr, DORMBR dormbr, XERBLA xerbla, DLAMCH dlamch, DLANGE dlange, ILAENV ilaenv, LSAME lsame)
        {


			#region Set Dependencies

			_dbdsdc = dbdsdc; _dgebrd = dgebrd; _dgelqf = dgelqf; _dgemm = dgemm; _dgeqrf = dgeqrf;
			_dlacpy = dlacpy; _dlascl = dlascl; _dlaset = dlaset; _dorgbr = dorgbr; _dorglq = dorglq;
			_dorgqr = dorgqr; _dormbr = dormbr; _xerbla = xerbla; _dlamch = dlamch; _dlange = dlange;
			_ilaenv = ilaenv; _lsame = lsame; 

            #endregion

        }
    
        public DGESDD()
        {
    

            #region Dependencies (Initialization)
            
            var lsame = new LSAME();
            var ieeeck = new IEEECK();
            var iparmq = new IPARMQ();
            var dlamc3 = new DLAMC3();
            var dlassq = new DLASSQ();
            var dcopy = new DCOPY();
            var xerbla = new XERBLA();
            var dlamrg = new DLAMRG();
            var dlapy2 = new DLAPY2();
            var drot = new DROT();
            var dnrm2 = new DNRM2();
            var dlasd5 = new DLASD5();
            var dlas2 = new DLAS2();
            var dlasq5 = new DLASQ5();
            var dlazq4 = new DLAZQ4();
            var dscal = new DSCAL();
            var dswap = new DSWAP();
            var dlasdt = new DLASDT();
            var ddot = new DDOT();
            var ilaenv = new ILAENV(ieeeck, iparmq);
            var dlamc1 = new DLAMC1(dlamc3);
            var dlamc4 = new DLAMC4(dlamc3);
            var dlamc5 = new DLAMC5(dlamc3);
            var dlamc2 = new DLAMC2(dlamc3, dlamc1, dlamc4, dlamc5);
            var dlamch = new DLAMCH(lsame, dlamc2);
            var dlanst = new DLANST(lsame, dlassq);
            var dlartg = new DLARTG(dlamch);
            var dlascl = new DLASCL(lsame, dlamch, xerbla);
            var dlacpy = new DLACPY(lsame);
            var dlaset = new DLASET(lsame);
            var dlasd2 = new DLASD2(dlamch, dlapy2, dcopy, dlacpy, dlamrg, dlaset, drot, xerbla);
            var dgemm = new DGEMM(lsame, xerbla);
            var dlaed6 = new DLAED6(dlamch);
            var dlasd4 = new DLASD4(dlaed6, dlasd5, dlamch);
            var dlasd3 = new DLASD3(dlamc3, dnrm2, dcopy, dgemm, dlacpy, dlascl, dlasd4, xerbla);
            var dlasd1 = new DLASD1(dlamrg, dlascl, dlasd2, dlasd3, xerbla);
            var dlasq6 = new DLASQ6(dlamch);
            var dlazq3 = new DLAZQ3(dlasq5, dlasq6, dlazq4, dlamch);
            var dlasrt = new DLASRT(lsame, xerbla);
            var dlasq2 = new DLASQ2(dlazq3, dlasrt, xerbla, dlamch, ilaenv);
            var dlasq1 = new DLASQ1(dcopy, dlas2, dlascl, dlasq2, dlasrt, xerbla, dlamch);
            var dlasr = new DLASR(lsame, xerbla);
            var dlasv2 = new DLASV2(dlamch);
            var dbdsqr = new DBDSQR(lsame, dlamch, dlartg, dlas2, dlasq1, dlasr, dlasv2, drot, dscal, dswap
                                       , xerbla);
            var dlasdq = new DLASDQ(dbdsqr, dlartg, dlasr, dswap, xerbla, lsame);
            var dlasd0 = new DLASD0(dlasd1, dlasdq, dlasdt, xerbla);
            var dlasd7 = new DLASD7(dcopy, dlamrg, drot, xerbla, dlamch, dlapy2);
            var dlasd8 = new DLASD8(dcopy, dlascl, dlasd4, dlaset, xerbla, ddot, dlamc3, dnrm2);
            var dlasd6 = new DLASD6(dcopy, dlamrg, dlascl, dlasd7, dlasd8, xerbla);
            var dlasda = new DLASDA(dcopy, dlasd6, dlasdq, dlasdt, dlaset, xerbla);
            var dbdsdc = new DBDSDC(lsame, ilaenv, dlamch, dlanst, dcopy, dlartg, dlascl, dlasd0, dlasda, dlasdq
                                       , dlaset, dlasr, dswap, xerbla);
            var dgemv = new DGEMV(lsame, xerbla);
            var dger = new DGER(xerbla);
            var dlarf = new DLARF(dgemv, dger, lsame);
            var dlarfg = new DLARFG(dlamch, dlapy2, dnrm2, dscal);
            var dgebd2 = new DGEBD2(dlarf, dlarfg, xerbla);
            var dlabrd = new DLABRD(dgemv, dlarfg, dscal);
            var dgebrd = new DGEBRD(dgebd2, dgemm, dlabrd, xerbla, ilaenv);
            var dgelq2 = new DGELQ2(dlarf, dlarfg, xerbla);
            var dtrmm = new DTRMM(lsame, xerbla);
            var dlarfb = new DLARFB(lsame, dcopy, dgemm, dtrmm);
            var dtrmv = new DTRMV(lsame, xerbla);
            var dlarft = new DLARFT(dgemv, dtrmv, lsame);
            var dgelqf = new DGELQF(dgelq2, dlarfb, dlarft, xerbla, ilaenv);
            var dgeqr2 = new DGEQR2(dlarf, dlarfg, xerbla);
            var dgeqrf = new DGEQRF(dgeqr2, dlarfb, dlarft, xerbla, ilaenv);
            var dorgl2 = new DORGL2(dlarf, dscal, xerbla);
            var dorglq = new DORGLQ(dlarfb, dlarft, dorgl2, xerbla, ilaenv);
            var dorg2r = new DORG2R(dlarf, dscal, xerbla);
            var dorgqr = new DORGQR(dlarfb, dlarft, dorg2r, xerbla, ilaenv);
            var dorgbr = new DORGBR(lsame, ilaenv, dorglq, dorgqr, xerbla);
            var dorml2 = new DORML2(lsame, dlarf, xerbla);
            var dormlq = new DORMLQ(lsame, ilaenv, dlarfb, dlarft, dorml2, xerbla);
            var dorm2r = new DORM2R(lsame, dlarf, xerbla);
            var dormqr = new DORMQR(lsame, ilaenv, dlarfb, dlarft, dorm2r, xerbla);
            var dormbr = new DORMBR(lsame, ilaenv, dormlq, dormqr, xerbla);
            var dlange = new DLANGE(dlassq, lsame);

			#endregion


			#region Set Dependencies

			_dbdsdc = dbdsdc; _dgebrd = dgebrd; _dgelqf = dgelqf; _dgemm = dgemm; _dgeqrf = dgeqrf;
			_dlacpy = dlacpy; _dlascl = dlascl; _dlaset = dlaset; _dorgbr = dorgbr; _dorglq = dorglq;
			_dorgqr = dorgqr; _dormbr = dormbr; _xerbla = xerbla; _dlamch = dlamch; _dlange = dlange;
			_ilaenv = ilaenv; _lsame = lsame; 

            #endregion

        }
        /// <summary>
        /// Purpose
        /// =======
        /// 
        /// DGESDD computes the singular value decomposition (SVD) of a real
        /// M-by-N matrix A, optionally computing the left and right singular
        /// vectors.  If singular vectors are desired, it uses a
        /// divide-and-conquer algorithm.
        /// 
        /// The SVD is written
        /// 
        /// A = U * SIGMA * transpose(V)
        /// 
        /// where SIGMA is an M-by-N matrix which is zero except for its
        /// min(m,n) diagonal elements, U is an M-by-M orthogonal matrix, and
        /// V is an N-by-N orthogonal matrix.  The diagonal elements of SIGMA
        /// are the singular values of A; they are real and non-negative, and
        /// are returned in descending order.  The first min(m,n) columns of
        /// U and V are the left and right singular vectors of A.
        /// 
        /// Note that the routine returns VT = V**T, not V.
        /// 
        /// The divide and conquer algorithm makes very mild assumptions about
        /// floating point arithmetic. It will work on machines with a guard
        /// digit in add/subtract, or on those binary machines without guard
        /// digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
        /// Cray-2. It could conceivably fail on hexadecimal or decimal machines
        /// without guard digits, but we know of none.
        /// 
        ///</summary>
        /// <param name="JOBZ">
        /// (input) CHARACTER*1
        /// Specifies options for computing all or part of the matrix U:
        /// = 'A':  all M columns of U and all N rows of V**T are
        /// returned in the arrays U and VT;
        /// = 'S':  the first min(M,N) columns of U and the first
        /// min(M,N) rows of V**T are returned in the arrays U
        /// and VT;
        /// = 'O':  If M .GE. N, the first N columns of U are overwritten
        /// on the array A and all rows of V**T are returned in
        /// the array VT;
        /// otherwise, all columns of U are returned in the
        /// array U and the first M rows of V**T are overwritten
        /// in the array A;
        /// = 'N':  no columns of U or rows of V**T are computed.
        ///</param>
        /// <param name="M">
        /// (input) INTEGER
        /// The number of rows of the input matrix A.  M .GE. 0.
        ///</param>
        /// <param name="N">
        /// (input) INTEGER
        /// The number of columns of the input matrix A.  N .GE. 0.
        ///</param>
        /// <param name="A">
        /// = U * SIGMA * transpose(V)
        ///</param>
        /// <param name="LDA">
        /// (input) INTEGER
        /// The leading dimension of the array A.  LDA .GE. max(1,M).
        ///</param>
        /// <param name="S">
        /// (output) DOUBLE PRECISION array, dimension (min(M,N))
        /// The singular values of A, sorted so that S(i) .GE. S(i+1).
        ///</param>
        /// <param name="U">
        /// (output) DOUBLE PRECISION array, dimension (LDU,UCOL)
        /// UCOL = M if JOBZ = 'A' or JOBZ = 'O' and M .LT. N;
        /// UCOL = min(M,N) if JOBZ = 'S'.
        /// If JOBZ = 'A' or JOBZ = 'O' and M .LT. N, U contains the M-by-M
        /// orthogonal matrix U;
        /// if JOBZ = 'S', U contains the first min(M,N) columns of U
        /// (the left singular vectors, stored columnwise);
        /// if JOBZ = 'O' and M .GE. N, or JOBZ = 'N', U is not referenced.
        ///</param>
        /// <param name="LDU">
        /// (input) INTEGER
        /// The leading dimension of the array U.  LDU .GE. 1; if
        /// JOBZ = 'S' or 'A' or JOBZ = 'O' and M .LT. N, LDU .GE. M.
        ///</param>
        /// <param name="VT">
        /// (output) DOUBLE PRECISION array, dimension (LDVT,N)
        /// If JOBZ = 'A' or JOBZ = 'O' and M .GE. N, VT contains the
        /// N-by-N orthogonal matrix V**T;
        /// if JOBZ = 'S', VT contains the first min(M,N) rows of
        /// V**T (the right singular vectors, stored rowwise);
        /// if JOBZ = 'O' and M .LT. N, or JOBZ = 'N', VT is not referenced.
        ///</param>
        /// <param name="LDVT">
        /// (input) INTEGER
        /// The leading dimension of the array VT.  LDVT .GE. 1; if
        /// JOBZ = 'A' or JOBZ = 'O' and M .GE. N, LDVT .GE. N;
        /// if JOBZ = 'S', LDVT .GE. min(M,N).
        ///</param>
        /// <param name="WORK">
        /// (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
        /// On exit, if INFO = 0, WORK(1) returns the optimal LWORK;
        ///</param>
        /// <param name="LWORK">
        /// (input) INTEGER
        /// The dimension of the array WORK. LWORK .GE. 1.
        /// If JOBZ = 'N',
        /// LWORK .GE. 3*min(M,N) + max(max(M,N),7*min(M,N)).
        /// If JOBZ = 'O',
        /// LWORK .GE. 3*min(M,N)*min(M,N) + 
        /// max(max(M,N),5*min(M,N)*min(M,N)+4*min(M,N)).
        /// If JOBZ = 'S' or 'A'
        /// LWORK .GE. 3*min(M,N)*min(M,N) +
        /// max(max(M,N),4*min(M,N)*min(M,N)+4*min(M,N)).
        /// For good performance, LWORK should generally be larger.
        /// If LWORK = -1 but other input arguments are legal, WORK(1)
        /// returns the optimal LWORK.
        ///</param>
        /// <param name="IWORK">
        /// (workspace) INTEGER array, dimension (8*min(M,N))
        ///</param>
        /// <param name="INFO">
        /// (output) INTEGER
        /// = 0:  successful exit.
        /// .LT. 0:  if INFO = -i, the i-th argument had an illegal value.
        /// .GT. 0:  DBDSDC did not converge, updating process failed.
        ///</param>
        public void Run(string JOBZ, int M, int N, ref double[] A, int offset_a, int LDA, ref double[] S, int offset_s
                         , ref double[] U, int offset_u, int LDU, ref double[] VT, int offset_vt, int LDVT, ref double[] WORK, int offset_work, int LWORK
                         , ref int[] IWORK, int offset_iwork, ref int INFO)
        {

            #region Variables
            
            var LQUERY = false; var WNTQA = false; var WNTQAS = false; var WNTQN = false; var WNTQO = false; 
            var WNTQS = false;var BDSPAC = 0; var BLK = 0; var CHUNK = 0; var I = 0; var IE = 0; var IERR = 0; var IL = 0; 
            var IR = 0;var ISCL = 0; var ITAU = 0; var ITAUP = 0; var ITAUQ = 0; var IU = 0; var IVT = 0; var LDWKVT = 0; 
            var LDWRKL = 0;var LDWRKR = 0; var LDWRKU = 0; var MAXWRK = 0; var MINMN = 0; var MINWRK = 0; var MNTHR = 0; 
            var NWORK = 0;var WRKBL = 0; double ANRM = 0; double BIGNUM = 0; double EPS = 0; double SMLNUM = 0; 
            var offset_idum = 0; var offset_dum = 0;

            #endregion


            #region Array Index Correction
            
             var o_a = -1 - LDA + offset_a;  var o_s = -1 + offset_s;  var o_u = -1 - LDU + offset_u; 
             var o_vt = -1 - LDVT + offset_vt; var o_work = -1 + offset_work;  var o_iwork = -1 + offset_iwork; 

            #endregion


            #region Strings
            
            JOBZ = JOBZ.Substring(0, 1);  

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
            // *  DGESDD computes the singular value decomposition (SVD) of a real
            // *  M-by-N matrix A, optionally computing the left and right singular
            // *  vectors.  If singular vectors are desired, it uses a
            // *  divide-and-conquer algorithm.
            // *
            // *  The SVD is written
            // *
            // *       A = U * SIGMA * transpose(V)
            // *
            // *  where SIGMA is an M-by-N matrix which is zero except for its
            // *  min(m,n) diagonal elements, U is an M-by-M orthogonal matrix, and
            // *  V is an N-by-N orthogonal matrix.  The diagonal elements of SIGMA
            // *  are the singular values of A; they are real and non-negative, and
            // *  are returned in descending order.  The first min(m,n) columns of
            // *  U and V are the left and right singular vectors of A.
            // *
            // *  Note that the routine returns VT = V**T, not V.
            // *
            // *  The divide and conquer algorithm makes very mild assumptions about
            // *  floating point arithmetic. It will work on machines with a guard
            // *  digit in add/subtract, or on those binary machines without guard
            // *  digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
            // *  Cray-2. It could conceivably fail on hexadecimal or decimal machines
            // *  without guard digits, but we know of none.
            // *
            // *  Arguments
            // *  =========
            // *
            // *  JOBZ    (input) CHARACTER*1
            // *          Specifies options for computing all or part of the matrix U:
            // *          = 'A':  all M columns of U and all N rows of V**T are
            // *                  returned in the arrays U and VT;
            // *          = 'S':  the first min(M,N) columns of U and the first
            // *                  min(M,N) rows of V**T are returned in the arrays U
            // *                  and VT;
            // *          = 'O':  If M >= N, the first N columns of U are overwritten
            // *                  on the array A and all rows of V**T are returned in
            // *                  the array VT;
            // *                  otherwise, all columns of U are returned in the
            // *                  array U and the first M rows of V**T are overwritten
            // *                  in the array A;
            // *          = 'N':  no columns of U or rows of V**T are computed.
            // *
            // *  M       (input) INTEGER
            // *          The number of rows of the input matrix A.  M >= 0.
            // *
            // *  N       (input) INTEGER
            // *          The number of columns of the input matrix A.  N >= 0.
            // *
            // *  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
            // *          On entry, the M-by-N matrix A.
            // *          On exit,
            // *          if JOBZ = 'O',  A is overwritten with the first N columns
            // *                          of U (the left singular vectors, stored
            // *                          columnwise) if M >= N;
            // *                          A is overwritten with the first M rows
            // *                          of V**T (the right singular vectors, stored
            // *                          rowwise) otherwise.
            // *          if JOBZ .ne. 'O', the contents of A are destroyed.
            // *
            // *  LDA     (input) INTEGER
            // *          The leading dimension of the array A.  LDA >= max(1,M).
            // *
            // *  S       (output) DOUBLE PRECISION array, dimension (min(M,N))
            // *          The singular values of A, sorted so that S(i) >= S(i+1).
            // *
            // *  U       (output) DOUBLE PRECISION array, dimension (LDU,UCOL)
            // *          UCOL = M if JOBZ = 'A' or JOBZ = 'O' and M < N;
            // *          UCOL = min(M,N) if JOBZ = 'S'.
            // *          If JOBZ = 'A' or JOBZ = 'O' and M < N, U contains the M-by-M
            // *          orthogonal matrix U;
            // *          if JOBZ = 'S', U contains the first min(M,N) columns of U
            // *          (the left singular vectors, stored columnwise);
            // *          if JOBZ = 'O' and M >= N, or JOBZ = 'N', U is not referenced.
            // *
            // *  LDU     (input) INTEGER
            // *          The leading dimension of the array U.  LDU >= 1; if
            // *          JOBZ = 'S' or 'A' or JOBZ = 'O' and M < N, LDU >= M.
            // *
            // *  VT      (output) DOUBLE PRECISION array, dimension (LDVT,N)
            // *          If JOBZ = 'A' or JOBZ = 'O' and M >= N, VT contains the
            // *          N-by-N orthogonal matrix V**T;
            // *          if JOBZ = 'S', VT contains the first min(M,N) rows of
            // *          V**T (the right singular vectors, stored rowwise);
            // *          if JOBZ = 'O' and M < N, or JOBZ = 'N', VT is not referenced.
            // *
            // *  LDVT    (input) INTEGER
            // *          The leading dimension of the array VT.  LDVT >= 1; if
            // *          JOBZ = 'A' or JOBZ = 'O' and M >= N, LDVT >= N;
            // *          if JOBZ = 'S', LDVT >= min(M,N).
            // *
            // *  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
            // *          On exit, if INFO = 0, WORK(1) returns the optimal LWORK;
            // *
            // *  LWORK   (input) INTEGER
            // *          The dimension of the array WORK. LWORK >= 1.
            // *          If JOBZ = 'N',
            // *            LWORK >= 3*min(M,N) + max(max(M,N),7*min(M,N)).
            // *          If JOBZ = 'O',
            // *            LWORK >= 3*min(M,N)*min(M,N) + 
            // *                     max(max(M,N),5*min(M,N)*min(M,N)+4*min(M,N)).
            // *          If JOBZ = 'S' or 'A'
            // *            LWORK >= 3*min(M,N)*min(M,N) +
            // *                     max(max(M,N),4*min(M,N)*min(M,N)+4*min(M,N)).
            // *          For good performance, LWORK should generally be larger.
            // *          If LWORK = -1 but other input arguments are legal, WORK(1)
            // *          returns the optimal LWORK.
            // *
            // *  IWORK   (workspace) INTEGER array, dimension (8*min(M,N))
            // *
            // *  INFO    (output) INTEGER
            // *          = 0:  successful exit.
            // *          < 0:  if INFO = -i, the i-th argument had an illegal value.
            // *          > 0:  DBDSDC did not converge, updating process failed.
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
            // *     .. Parameters ..
            // *     ..
            // *     .. Local Scalars ..
            // *     ..
            // *     .. Local Arrays ..
            // *     ..
            // *     .. External Subroutines ..
            // *     ..
            // *     .. External Functions ..
            // *     ..
            // *     .. Intrinsic Functions ..
            //      INTRINSIC          INT, MAX, MIN, SQRT;
            // *     ..
            // *     .. Executable Statements ..
            // *
            // *     Test the input arguments
            // *

            #endregion


            #region Body
            
            INFO = 0;
            MINMN = Math.Min(M, N);
            WNTQA = _lsame.Run(JOBZ, "A");
            WNTQS = _lsame.Run(JOBZ, "S");
            WNTQAS = WNTQA || WNTQS;
            WNTQO = _lsame.Run(JOBZ, "O");
            WNTQN = _lsame.Run(JOBZ, "N");
            LQUERY = LWORK ==  - 1;
            // *
            if (!(WNTQA || WNTQS || WNTQO || WNTQN))
            {
                INFO =  - 1;
            }
            else
            {
                if (M < 0)
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
                        if (LDA < Math.Max(1, M))
                        {
                            INFO =  - 5;
                        }
                        else
                        {
                            if (LDU < 1 || WNTQAS && LDU < M || WNTQO && M < N && LDU < M)
                            {
                                INFO =  - 8;
                            }
                            else
                            {
                                if (LDVT < 1 || WNTQA && LDVT < N || WNTQS && LDVT < MINMN || WNTQO && M >= N && LDVT < N)
                                {
                                    INFO =  - 10;
                                }
                            }
                        }
                    }
                }
            }
            // *
            // *     Compute workspace
            // *      (Note: Comments in the code beginning "Workspace:" describe the
            // *       minimal amount of workspace needed at that point in the code,
            // *       as well as the preferred amount for good performance.
            // *       NB refers to the optimal block size for the immediately
            // *       following subroutine, as returned by ILAENV.)
            // *
            if (INFO == 0)
            {
                MINWRK = 1;
                MAXWRK = 1;
                if (M >= N && MINMN > 0)
                {
                    // *
                    // *           Compute space needed for DBDSDC
                    // *
                    MNTHR = Convert.ToInt32(Math.Truncate(MINMN * 11.0E0 / 6.0E0));
                    if (WNTQN)
                    {
                        BDSPAC = 7 * N;
                    }
                    else
                    {
                        BDSPAC = 3 * N * N + 4 * N;
                    }
                    if (M >= MNTHR)
                    {
                        if (WNTQN)
                        {
                            // *
                            // *                 Path 1 (M much larger than N, JOBZ='N')
                            // *
                            WRKBL = N + N * _ilaenv.Run(1, "DGEQRF", " ", M, N,  - 1,  - 1);
                            WRKBL = Math.Max(WRKBL, 3 * N + 2 * N * _ilaenv.Run(1, "DGEBRD", " ", N, N,  - 1,  - 1));
                            MAXWRK = Math.Max(WRKBL, BDSPAC + N);
                            MINWRK = BDSPAC + N;
                        }
                        else
                        {
                            if (WNTQO)
                            {
                                // *
                                // *                 Path 2 (M much larger than N, JOBZ='O')
                                // *
                                WRKBL = N + N * _ilaenv.Run(1, "DGEQRF", " ", M, N,  - 1,  - 1);
                                WRKBL = Math.Max(WRKBL, N + N * _ilaenv.Run(1, "DORGQR", " ", M, N, N,  - 1));
                                WRKBL = Math.Max(WRKBL, 3 * N + 2 * N * _ilaenv.Run(1, "DGEBRD", " ", N, N,  - 1,  - 1));
                                WRKBL = Math.Max(WRKBL, 3 * N + N * _ilaenv.Run(1, "DORMBR", "QLN", N, N, N,  - 1));
                                WRKBL = Math.Max(WRKBL, 3 * N + N * _ilaenv.Run(1, "DORMBR", "PRT", N, N, N,  - 1));
                                WRKBL = Math.Max(WRKBL, BDSPAC + 3 * N);
                                MAXWRK = WRKBL + 2 * N * N;
                                MINWRK = BDSPAC + 2 * N * N + 3 * N;
                            }
                            else
                            {
                                if (WNTQS)
                                {
                                    // *
                                    // *                 Path 3 (M much larger than N, JOBZ='S')
                                    // *
                                    WRKBL = N + N * _ilaenv.Run(1, "DGEQRF", " ", M, N,  - 1,  - 1);
                                    WRKBL = Math.Max(WRKBL, N + N * _ilaenv.Run(1, "DORGQR", " ", M, N, N,  - 1));
                                    WRKBL = Math.Max(WRKBL, 3 * N + 2 * N * _ilaenv.Run(1, "DGEBRD", " ", N, N,  - 1,  - 1));
                                    WRKBL = Math.Max(WRKBL, 3 * N + N * _ilaenv.Run(1, "DORMBR", "QLN", N, N, N,  - 1));
                                    WRKBL = Math.Max(WRKBL, 3 * N + N * _ilaenv.Run(1, "DORMBR", "PRT", N, N, N,  - 1));
                                    WRKBL = Math.Max(WRKBL, BDSPAC + 3 * N);
                                    MAXWRK = WRKBL + N * N;
                                    MINWRK = BDSPAC + N * N + 3 * N;
                                }
                                else
                                {
                                    if (WNTQA)
                                    {
                                        // *
                                        // *                 Path 4 (M much larger than N, JOBZ='A')
                                        // *
                                        WRKBL = N + N * _ilaenv.Run(1, "DGEQRF", " ", M, N,  - 1,  - 1);
                                        WRKBL = Math.Max(WRKBL, N + M * _ilaenv.Run(1, "DORGQR", " ", M, M, N,  - 1));
                                        WRKBL = Math.Max(WRKBL, 3 * N + 2 * N * _ilaenv.Run(1, "DGEBRD", " ", N, N,  - 1,  - 1));
                                        WRKBL = Math.Max(WRKBL, 3 * N + N * _ilaenv.Run(1, "DORMBR", "QLN", N, N, N,  - 1));
                                        WRKBL = Math.Max(WRKBL, 3 * N + N * _ilaenv.Run(1, "DORMBR", "PRT", N, N, N,  - 1));
                                        WRKBL = Math.Max(WRKBL, BDSPAC + 3 * N);
                                        MAXWRK = WRKBL + N * N;
                                        MINWRK = BDSPAC + N * N + 3 * N;
                                    }
                                }
                            }
                        }
                    }
                    else
                    {
                        // *
                        // *              Path 5 (M at least N, but not much larger)
                        // *
                        WRKBL = 3 * N + (M + N) * _ilaenv.Run(1, "DGEBRD", " ", M, N,  - 1,  - 1);
                        if (WNTQN)
                        {
                            MAXWRK = Math.Max(WRKBL, BDSPAC + 3 * N);
                            MINWRK = 3 * N + Math.Max(M, BDSPAC);
                        }
                        else
                        {
                            if (WNTQO)
                            {
                                WRKBL = Math.Max(WRKBL, 3 * N + N * _ilaenv.Run(1, "DORMBR", "QLN", M, N, N,  - 1));
                                WRKBL = Math.Max(WRKBL, 3 * N + N * _ilaenv.Run(1, "DORMBR", "PRT", N, N, N,  - 1));
                                WRKBL = Math.Max(WRKBL, BDSPAC + 3 * N);
                                MAXWRK = WRKBL + M * N;
                                MINWRK = 3 * N + Math.Max(M, N * N + BDSPAC);
                            }
                            else
                            {
                                if (WNTQS)
                                {
                                    WRKBL = Math.Max(WRKBL, 3 * N + N * _ilaenv.Run(1, "DORMBR", "QLN", M, N, N,  - 1));
                                    WRKBL = Math.Max(WRKBL, 3 * N + N * _ilaenv.Run(1, "DORMBR", "PRT", N, N, N,  - 1));
                                    MAXWRK = Math.Max(WRKBL, BDSPAC + 3 * N);
                                    MINWRK = 3 * N + Math.Max(M, BDSPAC);
                                }
                                else
                                {
                                    if (WNTQA)
                                    {
                                        WRKBL = Math.Max(WRKBL, 3 * N + M * _ilaenv.Run(1, "DORMBR", "QLN", M, M, N,  - 1));
                                        WRKBL = Math.Max(WRKBL, 3 * N + N * _ilaenv.Run(1, "DORMBR", "PRT", N, N, N,  - 1));
                                        MAXWRK = Math.Max(MAXWRK, BDSPAC + 3 * N);
                                        MINWRK = 3 * N + Math.Max(M, BDSPAC);
                                    }
                                }
                            }
                        }
                    }
                }
                else
                {
                    if (MINMN > 0)
                    {
                        // *
                        // *           Compute space needed for DBDSDC
                        // *
                        MNTHR = Convert.ToInt32(Math.Truncate(MINMN * 11.0E0 / 6.0E0));
                        if (WNTQN)
                        {
                            BDSPAC = 7 * M;
                        }
                        else
                        {
                            BDSPAC = 3 * M * M + 4 * M;
                        }
                        if (N >= MNTHR)
                        {
                            if (WNTQN)
                            {
                                // *
                                // *                 Path 1t (N much larger than M, JOBZ='N')
                                // *
                                WRKBL = M + M * _ilaenv.Run(1, "DGELQF", " ", M, N,  - 1,  - 1);
                                WRKBL = Math.Max(WRKBL, 3 * M + 2 * M * _ilaenv.Run(1, "DGEBRD", " ", M, M,  - 1,  - 1));
                                MAXWRK = Math.Max(WRKBL, BDSPAC + M);
                                MINWRK = BDSPAC + M;
                            }
                            else
                            {
                                if (WNTQO)
                                {
                                    // *
                                    // *                 Path 2t (N much larger than M, JOBZ='O')
                                    // *
                                    WRKBL = M + M * _ilaenv.Run(1, "DGELQF", " ", M, N,  - 1,  - 1);
                                    WRKBL = Math.Max(WRKBL, M + M * _ilaenv.Run(1, "DORGLQ", " ", M, N, M,  - 1));
                                    WRKBL = Math.Max(WRKBL, 3 * M + 2 * M * _ilaenv.Run(1, "DGEBRD", " ", M, M,  - 1,  - 1));
                                    WRKBL = Math.Max(WRKBL, 3 * M + M * _ilaenv.Run(1, "DORMBR", "QLN", M, M, M,  - 1));
                                    WRKBL = Math.Max(WRKBL, 3 * M + M * _ilaenv.Run(1, "DORMBR", "PRT", M, M, M,  - 1));
                                    WRKBL = Math.Max(WRKBL, BDSPAC + 3 * M);
                                    MAXWRK = WRKBL + 2 * M * M;
                                    MINWRK = BDSPAC + 2 * M * M + 3 * M;
                                }
                                else
                                {
                                    if (WNTQS)
                                    {
                                        // *
                                        // *                 Path 3t (N much larger than M, JOBZ='S')
                                        // *
                                        WRKBL = M + M * _ilaenv.Run(1, "DGELQF", " ", M, N,  - 1,  - 1);
                                        WRKBL = Math.Max(WRKBL, M + M * _ilaenv.Run(1, "DORGLQ", " ", M, N, M,  - 1));
                                        WRKBL = Math.Max(WRKBL, 3 * M + 2 * M * _ilaenv.Run(1, "DGEBRD", " ", M, M,  - 1,  - 1));
                                        WRKBL = Math.Max(WRKBL, 3 * M + M * _ilaenv.Run(1, "DORMBR", "QLN", M, M, M,  - 1));
                                        WRKBL = Math.Max(WRKBL, 3 * M + M * _ilaenv.Run(1, "DORMBR", "PRT", M, M, M,  - 1));
                                        WRKBL = Math.Max(WRKBL, BDSPAC + 3 * M);
                                        MAXWRK = WRKBL + M * M;
                                        MINWRK = BDSPAC + M * M + 3 * M;
                                    }
                                    else
                                    {
                                        if (WNTQA)
                                        {
                                            // *
                                            // *                 Path 4t (N much larger than M, JOBZ='A')
                                            // *
                                            WRKBL = M + M * _ilaenv.Run(1, "DGELQF", " ", M, N,  - 1,  - 1);
                                            WRKBL = Math.Max(WRKBL, M + N * _ilaenv.Run(1, "DORGLQ", " ", N, N, M,  - 1));
                                            WRKBL = Math.Max(WRKBL, 3 * M + 2 * M * _ilaenv.Run(1, "DGEBRD", " ", M, M,  - 1,  - 1));
                                            WRKBL = Math.Max(WRKBL, 3 * M + M * _ilaenv.Run(1, "DORMBR", "QLN", M, M, M,  - 1));
                                            WRKBL = Math.Max(WRKBL, 3 * M + M * _ilaenv.Run(1, "DORMBR", "PRT", M, M, M,  - 1));
                                            WRKBL = Math.Max(WRKBL, BDSPAC + 3 * M);
                                            MAXWRK = WRKBL + M * M;
                                            MINWRK = BDSPAC + M * M + 3 * M;
                                        }
                                    }
                                }
                            }
                        }
                        else
                        {
                            // *
                            // *              Path 5t (N greater than M, but not much larger)
                            // *
                            WRKBL = 3 * M + (M + N) * _ilaenv.Run(1, "DGEBRD", " ", M, N,  - 1,  - 1);
                            if (WNTQN)
                            {
                                MAXWRK = Math.Max(WRKBL, BDSPAC + 3 * M);
                                MINWRK = 3 * M + Math.Max(N, BDSPAC);
                            }
                            else
                            {
                                if (WNTQO)
                                {
                                    WRKBL = Math.Max(WRKBL, 3 * M + M * _ilaenv.Run(1, "DORMBR", "QLN", M, M, N,  - 1));
                                    WRKBL = Math.Max(WRKBL, 3 * M + M * _ilaenv.Run(1, "DORMBR", "PRT", M, N, M,  - 1));
                                    WRKBL = Math.Max(WRKBL, BDSPAC + 3 * M);
                                    MAXWRK = WRKBL + M * N;
                                    MINWRK = 3 * M + Math.Max(N, M * M + BDSPAC);
                                }
                                else
                                {
                                    if (WNTQS)
                                    {
                                        WRKBL = Math.Max(WRKBL, 3 * M + M * _ilaenv.Run(1, "DORMBR", "QLN", M, M, N,  - 1));
                                        WRKBL = Math.Max(WRKBL, 3 * M + M * _ilaenv.Run(1, "DORMBR", "PRT", M, N, M,  - 1));
                                        MAXWRK = Math.Max(WRKBL, BDSPAC + 3 * M);
                                        MINWRK = 3 * M + Math.Max(N, BDSPAC);
                                    }
                                    else
                                    {
                                        if (WNTQA)
                                        {
                                            WRKBL = Math.Max(WRKBL, 3 * M + M * _ilaenv.Run(1, "DORMBR", "QLN", M, M, N,  - 1));
                                            WRKBL = Math.Max(WRKBL, 3 * M + M * _ilaenv.Run(1, "DORMBR", "PRT", N, N, M,  - 1));
                                            MAXWRK = Math.Max(WRKBL, BDSPAC + 3 * M);
                                            MINWRK = 3 * M + Math.Max(N, BDSPAC);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                MAXWRK = Math.Max(MAXWRK, MINWRK);
                WORK[1 + o_work] = MAXWRK;
                // *
                if (LWORK < MINWRK && !LQUERY)
                {
                    INFO =  - 12;
                }
            }
            // *
            if (INFO != 0)
            {
				_xerbla.Run("DGESDD",  - INFO);
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
            if (M == 0 || N == 0)
            {
                return;
            }
            // *
            // *     Get machine constants
            // *
            EPS = _dlamch.Run("P");
            SMLNUM = Math.Sqrt(_dlamch.Run("S")) / EPS;
            BIGNUM = ONE / SMLNUM;
            // *
            // *     Scale A if max element outside range [SMLNUM,BIGNUM]
            // *
            ANRM = _dlange.Run("M", M, N, A, offset_a, LDA, ref DUM, offset_dum);
            ISCL = 0;
            if (ANRM > ZERO && ANRM < SMLNUM)
            {
                ISCL = 1;
				_dlascl.Run("G", 0, 0, ANRM, SMLNUM, M
                                 , N, ref A, offset_a, LDA, ref IERR);
            }
            else
            {
                if (ANRM > BIGNUM)
                {
                    ISCL = 1;
					_dlascl.Run("G", 0, 0, ANRM, BIGNUM, M
                                     , N, ref A, offset_a, LDA, ref IERR);
                }
            }
            // *
            if (M >= N)
            {
                // *
                // *        A has at least as many rows as columns. If A has sufficiently
                // *        more rows than columns, first reduce using the QR
                // *        decomposition (if sufficient workspace available)
                // *
                if (M >= MNTHR)
                {
                    // *
                    if (WNTQN)
                    {
                        // *
                        // *              Path 1 (M much larger than N, JOBZ='N')
                        // *              No singular vectors to be computed
                        // *
                        ITAU = 1;
                        NWORK = ITAU + N;
						// *
						// *              Compute A=Q*R
						// *              (Workspace: need 2*N, prefer N+N*NB)
						// *
						_dgeqrf.Run(M, N, ref A, offset_a, LDA, ref WORK, ITAU + o_work, ref WORK, NWORK + o_work
                                         , LWORK - NWORK + 1, ref IERR);
						// *
						// *              Zero out below R
						// *
						_dlaset.Run("L", N - 1, N - 1, ZERO, ZERO, ref A, 2+1 * LDA + o_a
                                         , LDA);
                        IE = 1;
                        ITAUQ = IE + N;
                        ITAUP = ITAUQ + N;
                        NWORK = ITAUP + N;
						// *
						// *              Bidiagonalize R in A
						// *              (Workspace: need 4*N, prefer 3*N+2*N*NB)
						// *
						_dgebrd.Run(N, N, ref A, offset_a, LDA, ref S, offset_s, ref WORK, IE + o_work
                                         , ref WORK, ITAUQ + o_work, ref WORK, ITAUP + o_work, ref WORK, NWORK + o_work, LWORK - NWORK + 1, ref IERR);
                        NWORK = IE + N;
						// *
						// *              Perform bidiagonal SVD, computing singular values only
						// *              (Workspace: need N+BDSPAC)
						// *
						_dbdsdc.Run("U", "N", N, ref S, offset_s, ref WORK, IE + o_work, ref DUM, offset_dum
                                         , 1, ref DUM, offset_dum, 1, ref DUM, offset_dum, ref IDUM, offset_idum, ref WORK, NWORK + o_work
                                         , ref IWORK, offset_iwork, ref INFO);
                        // *
                    }
                    else
                    {
                        if (WNTQO)
                        {
                            // *
                            // *              Path 2 (M much larger than N, JOBZ = 'O')
                            // *              N left singular vectors to be overwritten on A and
                            // *              N right singular vectors to be computed in VT
                            // *
                            IR = 1;
                            // *
                            // *              WORK(IR) is LDWRKR by N
                            // *
                            if (LWORK >= LDA * N + N * N + 3 * N + BDSPAC)
                            {
                                LDWRKR = LDA;
                            }
                            else
                            {
                                LDWRKR = (LWORK - N * N - 3 * N - BDSPAC) / N;
                            }
                            ITAU = IR + LDWRKR * N;
                            NWORK = ITAU + N;
							// *
							// *              Compute A=Q*R
							// *              (Workspace: need N*N+2*N, prefer N*N+N+N*NB)
							// *
							_dgeqrf.Run(M, N, ref A, offset_a, LDA, ref WORK, ITAU + o_work, ref WORK, NWORK + o_work
                                             , LWORK - NWORK + 1, ref IERR);
							// *
							// *              Copy R to WORK(IR), zeroing out below it
							// *
							_dlacpy.Run("U", N, N, A, offset_a, LDA, ref WORK, IR + o_work
                                             , LDWRKR);
							_dlaset.Run("L", N - 1, N - 1, ZERO, ZERO, ref WORK, IR + 1 + o_work
                                             , LDWRKR);
							// *
							// *              Generate Q in A
							// *              (Workspace: need N*N+2*N, prefer N*N+N+N*NB)
							// *
							_dorgqr.Run(M, N, N, ref A, offset_a, LDA, WORK, ITAU + o_work
                                             , ref WORK, NWORK + o_work, LWORK - NWORK + 1, ref IERR);
                            IE = ITAU;
                            ITAUQ = IE + N;
                            ITAUP = ITAUQ + N;
                            NWORK = ITAUP + N;
							// *
							// *              Bidiagonalize R in VT, copying result to WORK(IR)
							// *              (Workspace: need N*N+4*N, prefer N*N+3*N+2*N*NB)
							// *
							_dgebrd.Run(N, N, ref WORK, IR + o_work, LDWRKR, ref S, offset_s, ref WORK, IE + o_work
                                             , ref WORK, ITAUQ + o_work, ref WORK, ITAUP + o_work, ref WORK, NWORK + o_work, LWORK - NWORK + 1, ref IERR);
                            // *
                            // *              WORK(IU) is N by N
                            // *
                            IU = NWORK;
                            NWORK = IU + N * N;
							// *
							// *              Perform bidiagonal SVD, computing left singular vectors
							// *              of bidiagonal matrix in WORK(IU) and computing right
							// *              singular vectors of bidiagonal matrix in VT
							// *              (Workspace: need N+N*N+BDSPAC)
							// *
							_dbdsdc.Run("U", "I", N, ref S, offset_s, ref WORK, IE + o_work, ref WORK, IU + o_work
                                             , N, ref VT, offset_vt, LDVT, ref DUM, offset_dum, ref IDUM, offset_idum, ref WORK, NWORK + o_work
                                             , ref IWORK, offset_iwork, ref INFO);
							// *
							// *              Overwrite WORK(IU) by left singular vectors of R
							// *              and VT by right singular vectors of R
							// *              (Workspace: need 2*N*N+3*N, prefer 2*N*N+2*N+N*NB)
							// *
							_dormbr.Run("Q", "L", "N", N, N, N
                                             , ref WORK, IR + o_work, LDWRKR, WORK, ITAUQ + o_work, ref WORK, IU + o_work, N, ref WORK, NWORK + o_work
                                             , LWORK - NWORK + 1, ref IERR);
							_dormbr.Run("P", "R", "T", N, N, N
                                             , ref WORK, IR + o_work, LDWRKR, WORK, ITAUP + o_work, ref VT, offset_vt, LDVT, ref WORK, NWORK + o_work
                                             , LWORK - NWORK + 1, ref IERR);
                            // *
                            // *              Multiply Q in A by left singular vectors of R in
                            // *              WORK(IU), storing result in WORK(IR) and copying to A
                            // *              (Workspace: need 2*N*N, prefer N*N+M*N)
                            // *
                            for (I = 1; LDWRKR >= 0 ? I <= M : I >= M; I += LDWRKR)
                            {
                                CHUNK = Math.Min(M - I + 1, LDWRKR);
								_dgemm.Run("N", "N", CHUNK, N, N, ONE
                                                , A, I+1 * LDA + o_a, LDA, WORK, IU + o_work, N, ZERO, ref WORK, IR + o_work
                                                , LDWRKR);
								_dlacpy.Run("F", CHUNK, N, WORK, IR + o_work, LDWRKR, ref A, I+1 * LDA + o_a
                                                 , LDA);
                            }
                            // *
                        }
                        else
                        {
                            if (WNTQS)
                            {
                                // *
                                // *              Path 3 (M much larger than N, JOBZ='S')
                                // *              N left singular vectors to be computed in U and
                                // *              N right singular vectors to be computed in VT
                                // *
                                IR = 1;
                                // *
                                // *              WORK(IR) is N by N
                                // *
                                LDWRKR = N;
                                ITAU = IR + LDWRKR * N;
                                NWORK = ITAU + N;
								// *
								// *              Compute A=Q*R
								// *              (Workspace: need N*N+2*N, prefer N*N+N+N*NB)
								// *
								_dgeqrf.Run(M, N, ref A, offset_a, LDA, ref WORK, ITAU + o_work, ref WORK, NWORK + o_work
                                                 , LWORK - NWORK + 1, ref IERR);
								// *
								// *              Copy R to WORK(IR), zeroing out below it
								// *
								_dlacpy.Run("U", N, N, A, offset_a, LDA, ref WORK, IR + o_work
                                                 , LDWRKR);
								_dlaset.Run("L", N - 1, N - 1, ZERO, ZERO, ref WORK, IR + 1 + o_work
                                                 , LDWRKR);
								// *
								// *              Generate Q in A
								// *              (Workspace: need N*N+2*N, prefer N*N+N+N*NB)
								// *
								_dorgqr.Run(M, N, N, ref A, offset_a, LDA, WORK, ITAU + o_work
                                                 , ref WORK, NWORK + o_work, LWORK - NWORK + 1, ref IERR);
                                IE = ITAU;
                                ITAUQ = IE + N;
                                ITAUP = ITAUQ + N;
                                NWORK = ITAUP + N;
								// *
								// *              Bidiagonalize R in WORK(IR)
								// *              (Workspace: need N*N+4*N, prefer N*N+3*N+2*N*NB)
								// *
								_dgebrd.Run(N, N, ref WORK, IR + o_work, LDWRKR, ref S, offset_s, ref WORK, IE + o_work
                                                 , ref WORK, ITAUQ + o_work, ref WORK, ITAUP + o_work, ref WORK, NWORK + o_work, LWORK - NWORK + 1, ref IERR);
								// *
								// *              Perform bidiagonal SVD, computing left singular vectors
								// *              of bidiagoal matrix in U and computing right singular
								// *              vectors of bidiagonal matrix in VT
								// *              (Workspace: need N+BDSPAC)
								// *
								_dbdsdc.Run("U", "I", N, ref S, offset_s, ref WORK, IE + o_work, ref U, offset_u
                                                 , LDU, ref VT, offset_vt, LDVT, ref DUM, offset_dum, ref IDUM, offset_idum, ref WORK, NWORK + o_work
                                                 , ref IWORK, offset_iwork, ref INFO);
								// *
								// *              Overwrite U by left singular vectors of R and VT
								// *              by right singular vectors of R
								// *              (Workspace: need N*N+3*N, prefer N*N+2*N+N*NB)
								// *
								_dormbr.Run("Q", "L", "N", N, N, N
                                                 , ref WORK, IR + o_work, LDWRKR, WORK, ITAUQ + o_work, ref U, offset_u, LDU, ref WORK, NWORK + o_work
                                                 , LWORK - NWORK + 1, ref IERR);
								// *
								_dormbr.Run("P", "R", "T", N, N, N
                                                 , ref WORK, IR + o_work, LDWRKR, WORK, ITAUP + o_work, ref VT, offset_vt, LDVT, ref WORK, NWORK + o_work
                                                 , LWORK - NWORK + 1, ref IERR);
								// *
								// *              Multiply Q in A by left singular vectors of R in
								// *              WORK(IR), storing result in U
								// *              (Workspace: need N*N)
								// *
								_dlacpy.Run("F", N, N, U, offset_u, LDU, ref WORK, IR + o_work
                                                 , LDWRKR);
								_dgemm.Run("N", "N", M, N, N, ONE
                                                , A, offset_a, LDA, WORK, IR + o_work, LDWRKR, ZERO, ref U, offset_u
                                                , LDU);
                                // *
                            }
                            else
                            {
                                if (WNTQA)
                                {
                                    // *
                                    // *              Path 4 (M much larger than N, JOBZ='A')
                                    // *              M left singular vectors to be computed in U and
                                    // *              N right singular vectors to be computed in VT
                                    // *
                                    IU = 1;
                                    // *
                                    // *              WORK(IU) is N by N
                                    // *
                                    LDWRKU = N;
                                    ITAU = IU + LDWRKU * N;
                                    NWORK = ITAU + N;
									// *
									// *              Compute A=Q*R, copying result to U
									// *              (Workspace: need N*N+2*N, prefer N*N+N+N*NB)
									// *
									_dgeqrf.Run(M, N, ref A, offset_a, LDA, ref WORK, ITAU + o_work, ref WORK, NWORK + o_work
                                                     , LWORK - NWORK + 1, ref IERR);
									_dlacpy.Run("L", M, N, A, offset_a, LDA, ref U, offset_u
                                                     , LDU);
									// *
									// *              Generate Q in U
									// *              (Workspace: need N*N+2*N, prefer N*N+N+N*NB)
									_dorgqr.Run(M, M, N, ref U, offset_u, LDU, WORK, ITAU + o_work
                                                     , ref WORK, NWORK + o_work, LWORK - NWORK + 1, ref IERR);
									// *
									// *              Produce R in A, zeroing out other entries
									// *
									_dlaset.Run("L", N - 1, N - 1, ZERO, ZERO, ref A, 2+1 * LDA + o_a
                                                     , LDA);
                                    IE = ITAU;
                                    ITAUQ = IE + N;
                                    ITAUP = ITAUQ + N;
                                    NWORK = ITAUP + N;
									// *
									// *              Bidiagonalize R in A
									// *              (Workspace: need N*N+4*N, prefer N*N+3*N+2*N*NB)
									// *
									_dgebrd.Run(N, N, ref A, offset_a, LDA, ref S, offset_s, ref WORK, IE + o_work
                                                     , ref WORK, ITAUQ + o_work, ref WORK, ITAUP + o_work, ref WORK, NWORK + o_work, LWORK - NWORK + 1, ref IERR);
									// *
									// *              Perform bidiagonal SVD, computing left singular vectors
									// *              of bidiagonal matrix in WORK(IU) and computing right
									// *              singular vectors of bidiagonal matrix in VT
									// *              (Workspace: need N+N*N+BDSPAC)
									// *
									_dbdsdc.Run("U", "I", N, ref S, offset_s, ref WORK, IE + o_work, ref WORK, IU + o_work
                                                     , N, ref VT, offset_vt, LDVT, ref DUM, offset_dum, ref IDUM, offset_idum, ref WORK, NWORK + o_work
                                                     , ref IWORK, offset_iwork, ref INFO);
									// *
									// *              Overwrite WORK(IU) by left singular vectors of R and VT
									// *              by right singular vectors of R
									// *              (Workspace: need N*N+3*N, prefer N*N+2*N+N*NB)
									// *
									_dormbr.Run("Q", "L", "N", N, N, N
                                                     , ref A, offset_a, LDA, WORK, ITAUQ + o_work, ref WORK, IU + o_work, LDWRKU, ref WORK, NWORK + o_work
                                                     , LWORK - NWORK + 1, ref IERR);
									_dormbr.Run("P", "R", "T", N, N, N
                                                     , ref A, offset_a, LDA, WORK, ITAUP + o_work, ref VT, offset_vt, LDVT, ref WORK, NWORK + o_work
                                                     , LWORK - NWORK + 1, ref IERR);
									// *
									// *              Multiply Q in U by left singular vectors of R in
									// *              WORK(IU), storing result in A
									// *              (Workspace: need N*N)
									// *
									_dgemm.Run("N", "N", M, N, N, ONE
                                                    , U, offset_u, LDU, WORK, IU + o_work, LDWRKU, ZERO, ref A, offset_a
                                                    , LDA);
									// *
									// *              Copy left singular vectors of A from A to U
									// *
									_dlacpy.Run("F", M, N, A, offset_a, LDA, ref U, offset_u
                                                     , LDU);
                                    // *
                                }
                            }
                        }
                    }
                    // *
                }
                else
                {
                    // *
                    // *           M .LT. MNTHR
                    // *
                    // *           Path 5 (M at least N, but not much larger)
                    // *           Reduce to bidiagonal form without QR decomposition
                    // *
                    IE = 1;
                    ITAUQ = IE + N;
                    ITAUP = ITAUQ + N;
                    NWORK = ITAUP + N;
					// *
					// *           Bidiagonalize A
					// *           (Workspace: need 3*N+M, prefer 3*N+(M+N)*NB)
					// *
					_dgebrd.Run(M, N, ref A, offset_a, LDA, ref S, offset_s, ref WORK, IE + o_work
                                     , ref WORK, ITAUQ + o_work, ref WORK, ITAUP + o_work, ref WORK, NWORK + o_work, LWORK - NWORK + 1, ref IERR);
                    if (WNTQN)
                    {
						// *
						// *              Perform bidiagonal SVD, only computing singular values
						// *              (Workspace: need N+BDSPAC)
						// *
						_dbdsdc.Run("U", "N", N, ref S, offset_s, ref WORK, IE + o_work, ref DUM, offset_dum
                                         , 1, ref DUM, offset_dum, 1, ref DUM, offset_dum, ref IDUM, offset_idum, ref WORK, NWORK + o_work
                                         , ref IWORK, offset_iwork, ref INFO);
                    }
                    else
                    {
                        if (WNTQO)
                        {
                            IU = NWORK;
                            if (LWORK >= M * N + 3 * N + BDSPAC)
                            {
                                // *
                                // *                 WORK( IU ) is M by N
                                // *
                                LDWRKU = M;
                                NWORK = IU + LDWRKU * N;
								_dlaset.Run("F", M, N, ZERO, ZERO, ref WORK, IU + o_work
                                                 , LDWRKU);
                            }
                            else
                            {
                                // *
                                // *                 WORK( IU ) is N by N
                                // *
                                LDWRKU = N;
                                NWORK = IU + LDWRKU * N;
                                // *
                                // *                 WORK(IR) is LDWRKR by N
                                // *
                                IR = NWORK;
                                LDWRKR = (LWORK - N * N - 3 * N) / N;
                            }
                            NWORK = IU + LDWRKU * N;
							// *
							// *              Perform bidiagonal SVD, computing left singular vectors
							// *              of bidiagonal matrix in WORK(IU) and computing right
							// *              singular vectors of bidiagonal matrix in VT
							// *              (Workspace: need N+N*N+BDSPAC)
							// *
							_dbdsdc.Run("U", "I", N, ref S, offset_s, ref WORK, IE + o_work, ref WORK, IU + o_work
                                             , LDWRKU, ref VT, offset_vt, LDVT, ref DUM, offset_dum, ref IDUM, offset_idum, ref WORK, NWORK + o_work
                                             , ref IWORK, offset_iwork, ref INFO);
							// *
							// *              Overwrite VT by right singular vectors of A
							// *              (Workspace: need N*N+2*N, prefer N*N+N+N*NB)
							// *
							_dormbr.Run("P", "R", "T", N, N, N
                                             , ref A, offset_a, LDA, WORK, ITAUP + o_work, ref VT, offset_vt, LDVT, ref WORK, NWORK + o_work
                                             , LWORK - NWORK + 1, ref IERR);
                            // *
                            if (LWORK >= M * N + 3 * N + BDSPAC)
                            {
								// *
								// *                 Overwrite WORK(IU) by left singular vectors of A
								// *                 (Workspace: need N*N+2*N, prefer N*N+N+N*NB)
								// *
								_dormbr.Run("Q", "L", "N", M, N, N
                                                 , ref A, offset_a, LDA, WORK, ITAUQ + o_work, ref WORK, IU + o_work, LDWRKU, ref WORK, NWORK + o_work
                                                 , LWORK - NWORK + 1, ref IERR);
								// *
								// *                 Copy left singular vectors of A from WORK(IU) to A
								// *
								_dlacpy.Run("F", M, N, WORK, IU + o_work, LDWRKU, ref A, offset_a
                                                 , LDA);
                            }
                            else
                            {
								// *
								// *                 Generate Q in A
								// *                 (Workspace: need N*N+2*N, prefer N*N+N+N*NB)
								// *
								_dorgbr.Run("Q", M, N, N, ref A, offset_a, LDA
                                                 , WORK, ITAUQ + o_work, ref WORK, NWORK + o_work, LWORK - NWORK + 1, ref IERR);
                                // *
                                // *                 Multiply Q in A by left singular vectors of
                                // *                 bidiagonal matrix in WORK(IU), storing result in
                                // *                 WORK(IR) and copying to A
                                // *                 (Workspace: need 2*N*N, prefer N*N+M*N)
                                // *
                                for (I = 1; LDWRKR >= 0 ? I <= M : I >= M; I += LDWRKR)
                                {
                                    CHUNK = Math.Min(M - I + 1, LDWRKR);
									_dgemm.Run("N", "N", CHUNK, N, N, ONE
                                                    , A, I+1 * LDA + o_a, LDA, WORK, IU + o_work, LDWRKU, ZERO, ref WORK, IR + o_work
                                                    , LDWRKR);
									_dlacpy.Run("F", CHUNK, N, WORK, IR + o_work, LDWRKR, ref A, I+1 * LDA + o_a
                                                     , LDA);
                                }
                            }
                            // *
                        }
                        else
                        {
                            if (WNTQS)
                            {
								// *
								// *              Perform bidiagonal SVD, computing left singular vectors
								// *              of bidiagonal matrix in U and computing right singular
								// *              vectors of bidiagonal matrix in VT
								// *              (Workspace: need N+BDSPAC)
								// *
								_dlaset.Run("F", M, N, ZERO, ZERO, ref U, offset_u
                                                 , LDU);
								_dbdsdc.Run("U", "I", N, ref S, offset_s, ref WORK, IE + o_work, ref U, offset_u
                                                 , LDU, ref VT, offset_vt, LDVT, ref DUM, offset_dum, ref IDUM, offset_idum, ref WORK, NWORK + o_work
                                                 , ref IWORK, offset_iwork, ref INFO);
								// *
								// *              Overwrite U by left singular vectors of A and VT
								// *              by right singular vectors of A
								// *              (Workspace: need 3*N, prefer 2*N+N*NB)
								// *
								_dormbr.Run("Q", "L", "N", M, N, N
                                                 , ref A, offset_a, LDA, WORK, ITAUQ + o_work, ref U, offset_u, LDU, ref WORK, NWORK + o_work
                                                 , LWORK - NWORK + 1, ref IERR);
								_dormbr.Run("P", "R", "T", N, N, N
                                                 , ref A, offset_a, LDA, WORK, ITAUP + o_work, ref VT, offset_vt, LDVT, ref WORK, NWORK + o_work
                                                 , LWORK - NWORK + 1, ref IERR);
                            }
                            else
                            {
                                if (WNTQA)
                                {
									// *
									// *              Perform bidiagonal SVD, computing left singular vectors
									// *              of bidiagonal matrix in U and computing right singular
									// *              vectors of bidiagonal matrix in VT
									// *              (Workspace: need N+BDSPAC)
									// *
									_dlaset.Run("F", M, M, ZERO, ZERO, ref U, offset_u
                                                     , LDU);
									_dbdsdc.Run("U", "I", N, ref S, offset_s, ref WORK, IE + o_work, ref U, offset_u
                                                     , LDU, ref VT, offset_vt, LDVT, ref DUM, offset_dum, ref IDUM, offset_idum, ref WORK, NWORK + o_work
                                                     , ref IWORK, offset_iwork, ref INFO);
                                    // *
                                    // *              Set the right corner of U to identity matrix
                                    // *
                                    if (M > N)
                                    {
										_dlaset.Run("F", M - N, M - N, ZERO, ONE, ref U, N + 1+(N + 1) * LDU + o_u
                                                         , LDU);
                                    }
									// *
									// *              Overwrite U by left singular vectors of A and VT
									// *              by right singular vectors of A
									// *              (Workspace: need N*N+2*N+M, prefer N*N+2*N+M*NB)
									// *
									_dormbr.Run("Q", "L", "N", M, M, N
                                                     , ref A, offset_a, LDA, WORK, ITAUQ + o_work, ref U, offset_u, LDU, ref WORK, NWORK + o_work
                                                     , LWORK - NWORK + 1, ref IERR);
									_dormbr.Run("P", "R", "T", N, N, M
                                                     , ref A, offset_a, LDA, WORK, ITAUP + o_work, ref VT, offset_vt, LDVT, ref WORK, NWORK + o_work
                                                     , LWORK - NWORK + 1, ref IERR);
                                }
                            }
                        }
                    }
                    // *
                }
                // *
            }
            else
            {
                // *
                // *        A has more columns than rows. If A has sufficiently more
                // *        columns than rows, first reduce using the LQ decomposition (if
                // *        sufficient workspace available)
                // *
                if (N >= MNTHR)
                {
                    // *
                    if (WNTQN)
                    {
                        // *
                        // *              Path 1t (N much larger than M, JOBZ='N')
                        // *              No singular vectors to be computed
                        // *
                        ITAU = 1;
                        NWORK = ITAU + M;
						// *
						// *              Compute A=L*Q
						// *              (Workspace: need 2*M, prefer M+M*NB)
						// *
						_dgelqf.Run(M, N, ref A, offset_a, LDA, ref WORK, ITAU + o_work, ref WORK, NWORK + o_work
                                         , LWORK - NWORK + 1, ref IERR);
						// *
						// *              Zero out above L
						// *
						_dlaset.Run("U", M - 1, M - 1, ZERO, ZERO, ref A, 1+2 * LDA + o_a
                                         , LDA);
                        IE = 1;
                        ITAUQ = IE + M;
                        ITAUP = ITAUQ + M;
                        NWORK = ITAUP + M;
						// *
						// *              Bidiagonalize L in A
						// *              (Workspace: need 4*M, prefer 3*M+2*M*NB)
						// *
						_dgebrd.Run(M, M, ref A, offset_a, LDA, ref S, offset_s, ref WORK, IE + o_work
                                         , ref WORK, ITAUQ + o_work, ref WORK, ITAUP + o_work, ref WORK, NWORK + o_work, LWORK - NWORK + 1, ref IERR);
                        NWORK = IE + M;
						// *
						// *              Perform bidiagonal SVD, computing singular values only
						// *              (Workspace: need M+BDSPAC)
						// *
						_dbdsdc.Run("U", "N", M, ref S, offset_s, ref WORK, IE + o_work, ref DUM, offset_dum
                                         , 1, ref DUM, offset_dum, 1, ref DUM, offset_dum, ref IDUM, offset_idum, ref WORK, NWORK + o_work
                                         , ref IWORK, offset_iwork, ref INFO);
                        // *
                    }
                    else
                    {
                        if (WNTQO)
                        {
                            // *
                            // *              Path 2t (N much larger than M, JOBZ='O')
                            // *              M right singular vectors to be overwritten on A and
                            // *              M left singular vectors to be computed in U
                            // *
                            IVT = 1;
                            // *
                            // *              IVT is M by M
                            // *
                            IL = IVT + M * M;
                            if (LWORK >= M * N + M * M + 3 * M + BDSPAC)
                            {
                                // *
                                // *                 WORK(IL) is M by N
                                // *
                                LDWRKL = M;
                                CHUNK = N;
                            }
                            else
                            {
                                LDWRKL = M;
                                CHUNK = (LWORK - M * M) / M;
                            }
                            ITAU = IL + LDWRKL * M;
                            NWORK = ITAU + M;
							// *
							// *              Compute A=L*Q
							// *              (Workspace: need M*M+2*M, prefer M*M+M+M*NB)
							// *
							_dgelqf.Run(M, N, ref A, offset_a, LDA, ref WORK, ITAU + o_work, ref WORK, NWORK + o_work
                                             , LWORK - NWORK + 1, ref IERR);
							// *
							// *              Copy L to WORK(IL), zeroing about above it
							// *
							_dlacpy.Run("L", M, M, A, offset_a, LDA, ref WORK, IL + o_work
                                             , LDWRKL);
							_dlaset.Run("U", M - 1, M - 1, ZERO, ZERO, ref WORK, IL + LDWRKL + o_work
                                             , LDWRKL);
							// *
							// *              Generate Q in A
							// *              (Workspace: need M*M+2*M, prefer M*M+M+M*NB)
							// *
							_dorglq.Run(M, N, M, ref A, offset_a, LDA, WORK, ITAU + o_work
                                             , ref WORK, NWORK + o_work, LWORK - NWORK + 1, ref IERR);
                            IE = ITAU;
                            ITAUQ = IE + M;
                            ITAUP = ITAUQ + M;
                            NWORK = ITAUP + M;
							// *
							// *              Bidiagonalize L in WORK(IL)
							// *              (Workspace: need M*M+4*M, prefer M*M+3*M+2*M*NB)
							// *
							_dgebrd.Run(M, M, ref WORK, IL + o_work, LDWRKL, ref S, offset_s, ref WORK, IE + o_work
                                             , ref WORK, ITAUQ + o_work, ref WORK, ITAUP + o_work, ref WORK, NWORK + o_work, LWORK - NWORK + 1, ref IERR);
							// *
							// *              Perform bidiagonal SVD, computing left singular vectors
							// *              of bidiagonal matrix in U, and computing right singular
							// *              vectors of bidiagonal matrix in WORK(IVT)
							// *              (Workspace: need M+M*M+BDSPAC)
							// *
							_dbdsdc.Run("U", "I", M, ref S, offset_s, ref WORK, IE + o_work, ref U, offset_u
                                             , LDU, ref WORK, IVT + o_work, M, ref DUM, offset_dum, ref IDUM, offset_idum, ref WORK, NWORK + o_work
                                             , ref IWORK, offset_iwork, ref INFO);
							// *
							// *              Overwrite U by left singular vectors of L and WORK(IVT)
							// *              by right singular vectors of L
							// *              (Workspace: need 2*M*M+3*M, prefer 2*M*M+2*M+M*NB)
							// *
							_dormbr.Run("Q", "L", "N", M, M, M
                                             , ref WORK, IL + o_work, LDWRKL, WORK, ITAUQ + o_work, ref U, offset_u, LDU, ref WORK, NWORK + o_work
                                             , LWORK - NWORK + 1, ref IERR);
							_dormbr.Run("P", "R", "T", M, M, M
                                             , ref WORK, IL + o_work, LDWRKL, WORK, ITAUP + o_work, ref WORK, IVT + o_work, M, ref WORK, NWORK + o_work
                                             , LWORK - NWORK + 1, ref IERR);
                            // *
                            // *              Multiply right singular vectors of L in WORK(IVT) by Q
                            // *              in A, storing result in WORK(IL) and copying to A
                            // *              (Workspace: need 2*M*M, prefer M*M+M*N)
                            // *
                            for (I = 1; CHUNK >= 0 ? I <= N : I >= N; I += CHUNK)
                            {
                                BLK = Math.Min(N - I + 1, CHUNK);
								_dgemm.Run("N", "N", M, BLK, M, ONE
                                                , WORK, IVT + o_work, M, A, 1+I * LDA + o_a, LDA, ZERO, ref WORK, IL + o_work
                                                , LDWRKL);
								_dlacpy.Run("F", M, BLK, WORK, IL + o_work, LDWRKL, ref A, 1+I * LDA + o_a
                                                 , LDA);
                            }
                            // *
                        }
                        else
                        {
                            if (WNTQS)
                            {
                                // *
                                // *              Path 3t (N much larger than M, JOBZ='S')
                                // *              M right singular vectors to be computed in VT and
                                // *              M left singular vectors to be computed in U
                                // *
                                IL = 1;
                                // *
                                // *              WORK(IL) is M by M
                                // *
                                LDWRKL = M;
                                ITAU = IL + LDWRKL * M;
                                NWORK = ITAU + M;
								// *
								// *              Compute A=L*Q
								// *              (Workspace: need M*M+2*M, prefer M*M+M+M*NB)
								// *
								_dgelqf.Run(M, N, ref A, offset_a, LDA, ref WORK, ITAU + o_work, ref WORK, NWORK + o_work
                                                 , LWORK - NWORK + 1, ref IERR);
								// *
								// *              Copy L to WORK(IL), zeroing out above it
								// *
								_dlacpy.Run("L", M, M, A, offset_a, LDA, ref WORK, IL + o_work
                                                 , LDWRKL);
								_dlaset.Run("U", M - 1, M - 1, ZERO, ZERO, ref WORK, IL + LDWRKL + o_work
                                                 , LDWRKL);
								// *
								// *              Generate Q in A
								// *              (Workspace: need M*M+2*M, prefer M*M+M+M*NB)
								// *
								_dorglq.Run(M, N, M, ref A, offset_a, LDA, WORK, ITAU + o_work
                                                 , ref WORK, NWORK + o_work, LWORK - NWORK + 1, ref IERR);
                                IE = ITAU;
                                ITAUQ = IE + M;
                                ITAUP = ITAUQ + M;
                                NWORK = ITAUP + M;
								// *
								// *              Bidiagonalize L in WORK(IU), copying result to U
								// *              (Workspace: need M*M+4*M, prefer M*M+3*M+2*M*NB)
								// *
								_dgebrd.Run(M, M, ref WORK, IL + o_work, LDWRKL, ref S, offset_s, ref WORK, IE + o_work
                                                 , ref WORK, ITAUQ + o_work, ref WORK, ITAUP + o_work, ref WORK, NWORK + o_work, LWORK - NWORK + 1, ref IERR);
								// *
								// *              Perform bidiagonal SVD, computing left singular vectors
								// *              of bidiagonal matrix in U and computing right singular
								// *              vectors of bidiagonal matrix in VT
								// *              (Workspace: need M+BDSPAC)
								// *
								_dbdsdc.Run("U", "I", M, ref S, offset_s, ref WORK, IE + o_work, ref U, offset_u
                                                 , LDU, ref VT, offset_vt, LDVT, ref DUM, offset_dum, ref IDUM, offset_idum, ref WORK, NWORK + o_work
                                                 , ref IWORK, offset_iwork, ref INFO);
								// *
								// *              Overwrite U by left singular vectors of L and VT
								// *              by right singular vectors of L
								// *              (Workspace: need M*M+3*M, prefer M*M+2*M+M*NB)
								// *
								_dormbr.Run("Q", "L", "N", M, M, M
                                                 , ref WORK, IL + o_work, LDWRKL, WORK, ITAUQ + o_work, ref U, offset_u, LDU, ref WORK, NWORK + o_work
                                                 , LWORK - NWORK + 1, ref IERR);
								_dormbr.Run("P", "R", "T", M, M, M
                                                 , ref WORK, IL + o_work, LDWRKL, WORK, ITAUP + o_work, ref VT, offset_vt, LDVT, ref WORK, NWORK + o_work
                                                 , LWORK - NWORK + 1, ref IERR);
								// *
								// *              Multiply right singular vectors of L in WORK(IL) by
								// *              Q in A, storing result in VT
								// *              (Workspace: need M*M)
								// *
								_dlacpy.Run("F", M, M, VT, offset_vt, LDVT, ref WORK, IL + o_work
                                                 , LDWRKL);
								_dgemm.Run("N", "N", M, N, M, ONE
                                                , WORK, IL + o_work, LDWRKL, A, offset_a, LDA, ZERO, ref VT, offset_vt
                                                , LDVT);
                                // *
                            }
                            else
                            {
                                if (WNTQA)
                                {
                                    // *
                                    // *              Path 4t (N much larger than M, JOBZ='A')
                                    // *              N right singular vectors to be computed in VT and
                                    // *              M left singular vectors to be computed in U
                                    // *
                                    IVT = 1;
                                    // *
                                    // *              WORK(IVT) is M by M
                                    // *
                                    LDWKVT = M;
                                    ITAU = IVT + LDWKVT * M;
                                    NWORK = ITAU + M;
									// *
									// *              Compute A=L*Q, copying result to VT
									// *              (Workspace: need M*M+2*M, prefer M*M+M+M*NB)
									// *
									_dgelqf.Run(M, N, ref A, offset_a, LDA, ref WORK, ITAU + o_work, ref WORK, NWORK + o_work
                                                     , LWORK - NWORK + 1, ref IERR);
									_dlacpy.Run("U", M, N, A, offset_a, LDA, ref VT, offset_vt
                                                     , LDVT);
									// *
									// *              Generate Q in VT
									// *              (Workspace: need M*M+2*M, prefer M*M+M+M*NB)
									// *
									_dorglq.Run(N, N, M, ref VT, offset_vt, LDVT, WORK, ITAU + o_work
                                                     , ref WORK, NWORK + o_work, LWORK - NWORK + 1, ref IERR);
									// *
									// *              Produce L in A, zeroing out other entries
									// *
									_dlaset.Run("U", M - 1, M - 1, ZERO, ZERO, ref A, 1+2 * LDA + o_a
                                                     , LDA);
                                    IE = ITAU;
                                    ITAUQ = IE + M;
                                    ITAUP = ITAUQ + M;
                                    NWORK = ITAUP + M;
									// *
									// *              Bidiagonalize L in A
									// *              (Workspace: need M*M+4*M, prefer M*M+3*M+2*M*NB)
									// *
									_dgebrd.Run(M, M, ref A, offset_a, LDA, ref S, offset_s, ref WORK, IE + o_work
                                                     , ref WORK, ITAUQ + o_work, ref WORK, ITAUP + o_work, ref WORK, NWORK + o_work, LWORK - NWORK + 1, ref IERR);
									// *
									// *              Perform bidiagonal SVD, computing left singular vectors
									// *              of bidiagonal matrix in U and computing right singular
									// *              vectors of bidiagonal matrix in WORK(IVT)
									// *              (Workspace: need M+M*M+BDSPAC)
									// *
									_dbdsdc.Run("U", "I", M, ref S, offset_s, ref WORK, IE + o_work, ref U, offset_u
                                                     , LDU, ref WORK, IVT + o_work, LDWKVT, ref DUM, offset_dum, ref IDUM, offset_idum, ref WORK, NWORK + o_work
                                                     , ref IWORK, offset_iwork, ref INFO);
									// *
									// *              Overwrite U by left singular vectors of L and WORK(IVT)
									// *              by right singular vectors of L
									// *              (Workspace: need M*M+3*M, prefer M*M+2*M+M*NB)
									// *
									_dormbr.Run("Q", "L", "N", M, M, M
                                                     , ref A, offset_a, LDA, WORK, ITAUQ + o_work, ref U, offset_u, LDU, ref WORK, NWORK + o_work
                                                     , LWORK - NWORK + 1, ref IERR);
									_dormbr.Run("P", "R", "T", M, M, M
                                                     , ref A, offset_a, LDA, WORK, ITAUP + o_work, ref WORK, IVT + o_work, LDWKVT, ref WORK, NWORK + o_work
                                                     , LWORK - NWORK + 1, ref IERR);
									// *
									// *              Multiply right singular vectors of L in WORK(IVT) by
									// *              Q in VT, storing result in A
									// *              (Workspace: need M*M)
									// *
									_dgemm.Run("N", "N", M, N, M, ONE
                                                    , WORK, IVT + o_work, LDWKVT, VT, offset_vt, LDVT, ZERO, ref A, offset_a
                                                    , LDA);
									// *
									// *              Copy right singular vectors of A from A to VT
									// *
									_dlacpy.Run("F", M, N, A, offset_a, LDA, ref VT, offset_vt
                                                     , LDVT);
                                    // *
                                }
                            }
                        }
                    }
                    // *
                }
                else
                {
                    // *
                    // *           N .LT. MNTHR
                    // *
                    // *           Path 5t (N greater than M, but not much larger)
                    // *           Reduce to bidiagonal form without LQ decomposition
                    // *
                    IE = 1;
                    ITAUQ = IE + M;
                    ITAUP = ITAUQ + M;
                    NWORK = ITAUP + M;
					// *
					// *           Bidiagonalize A
					// *           (Workspace: need 3*M+N, prefer 3*M+(M+N)*NB)
					// *
					_dgebrd.Run(M, N, ref A, offset_a, LDA, ref S, offset_s, ref WORK, IE + o_work
                                     , ref WORK, ITAUQ + o_work, ref WORK, ITAUP + o_work, ref WORK, NWORK + o_work, LWORK - NWORK + 1, ref IERR);
                    if (WNTQN)
                    {
						// *
						// *              Perform bidiagonal SVD, only computing singular values
						// *              (Workspace: need M+BDSPAC)
						// *
						_dbdsdc.Run("L", "N", M, ref S, offset_s, ref WORK, IE + o_work, ref DUM, offset_dum
                                         , 1, ref DUM, offset_dum, 1, ref DUM, offset_dum, ref IDUM, offset_idum, ref WORK, NWORK + o_work
                                         , ref IWORK, offset_iwork, ref INFO);
                    }
                    else
                    {
                        if (WNTQO)
                        {
                            LDWKVT = M;
                            IVT = NWORK;
                            if (LWORK >= M * N + 3 * M + BDSPAC)
                            {
								// *
								// *                 WORK( IVT ) is M by N
								// *
								_dlaset.Run("F", M, N, ZERO, ZERO, ref WORK, IVT + o_work
                                                 , LDWKVT);
                                NWORK = IVT + LDWKVT * N;
                            }
                            else
                            {
                                // *
                                // *                 WORK( IVT ) is M by M
                                // *
                                NWORK = IVT + LDWKVT * M;
                                IL = NWORK;
                                // *
                                // *                 WORK(IL) is M by CHUNK
                                // *
                                CHUNK = (LWORK - M * M - 3 * M) / M;
                            }
							// *
							// *              Perform bidiagonal SVD, computing left singular vectors
							// *              of bidiagonal matrix in U and computing right singular
							// *              vectors of bidiagonal matrix in WORK(IVT)
							// *              (Workspace: need M*M+BDSPAC)
							// *
							_dbdsdc.Run("L", "I", M, ref S, offset_s, ref WORK, IE + o_work, ref U, offset_u
                                             , LDU, ref WORK, IVT + o_work, LDWKVT, ref DUM, offset_dum, ref IDUM, offset_idum, ref WORK, NWORK + o_work
                                             , ref IWORK, offset_iwork, ref INFO);
							// *
							// *              Overwrite U by left singular vectors of A
							// *              (Workspace: need M*M+2*M, prefer M*M+M+M*NB)
							// *
							_dormbr.Run("Q", "L", "N", M, M, N
                                             , ref A, offset_a, LDA, WORK, ITAUQ + o_work, ref U, offset_u, LDU, ref WORK, NWORK + o_work
                                             , LWORK - NWORK + 1, ref IERR);
                            // *
                            if (LWORK >= M * N + 3 * M + BDSPAC)
                            {
								// *
								// *                 Overwrite WORK(IVT) by left singular vectors of A
								// *                 (Workspace: need M*M+2*M, prefer M*M+M+M*NB)
								// *
								_dormbr.Run("P", "R", "T", M, N, M
                                                 , ref A, offset_a, LDA, WORK, ITAUP + o_work, ref WORK, IVT + o_work, LDWKVT, ref WORK, NWORK + o_work
                                                 , LWORK - NWORK + 1, ref IERR);
								// *
								// *                 Copy right singular vectors of A from WORK(IVT) to A
								// *
								_dlacpy.Run("F", M, N, WORK, IVT + o_work, LDWKVT, ref A, offset_a
                                                 , LDA);
                            }
                            else
                            {
								// *
								// *                 Generate P**T in A
								// *                 (Workspace: need M*M+2*M, prefer M*M+M+M*NB)
								// *
								_dorgbr.Run("P", M, N, M, ref A, offset_a, LDA
                                                 , WORK, ITAUP + o_work, ref WORK, NWORK + o_work, LWORK - NWORK + 1, ref IERR);
                                // *
                                // *                 Multiply Q in A by right singular vectors of
                                // *                 bidiagonal matrix in WORK(IVT), storing result in
                                // *                 WORK(IL) and copying to A
                                // *                 (Workspace: need 2*M*M, prefer M*M+M*N)
                                // *
                                for (I = 1; CHUNK >= 0 ? I <= N : I >= N; I += CHUNK)
                                {
                                    BLK = Math.Min(N - I + 1, CHUNK);
									_dgemm.Run("N", "N", M, BLK, M, ONE
                                                    , WORK, IVT + o_work, LDWKVT, A, 1+I * LDA + o_a, LDA, ZERO, ref WORK, IL + o_work
                                                    , M);
									_dlacpy.Run("F", M, BLK, WORK, IL + o_work, M, ref A, 1+I * LDA + o_a
                                                     , LDA);
                                }
                            }
                        }
                        else
                        {
                            if (WNTQS)
                            {
								// *
								// *              Perform bidiagonal SVD, computing left singular vectors
								// *              of bidiagonal matrix in U and computing right singular
								// *              vectors of bidiagonal matrix in VT
								// *              (Workspace: need M+BDSPAC)
								// *
								_dlaset.Run("F", M, N, ZERO, ZERO, ref VT, offset_vt
                                                 , LDVT);
								_dbdsdc.Run("L", "I", M, ref S, offset_s, ref WORK, IE + o_work, ref U, offset_u
                                                 , LDU, ref VT, offset_vt, LDVT, ref DUM, offset_dum, ref IDUM, offset_idum, ref WORK, NWORK + o_work
                                                 , ref IWORK, offset_iwork, ref INFO);
								// *
								// *              Overwrite U by left singular vectors of A and VT
								// *              by right singular vectors of A
								// *              (Workspace: need 3*M, prefer 2*M+M*NB)
								// *
								_dormbr.Run("Q", "L", "N", M, M, N
                                                 , ref A, offset_a, LDA, WORK, ITAUQ + o_work, ref U, offset_u, LDU, ref WORK, NWORK + o_work
                                                 , LWORK - NWORK + 1, ref IERR);
								_dormbr.Run("P", "R", "T", M, N, M
                                                 , ref A, offset_a, LDA, WORK, ITAUP + o_work, ref VT, offset_vt, LDVT, ref WORK, NWORK + o_work
                                                 , LWORK - NWORK + 1, ref IERR);
                            }
                            else
                            {
                                if (WNTQA)
                                {
									// *
									// *              Perform bidiagonal SVD, computing left singular vectors
									// *              of bidiagonal matrix in U and computing right singular
									// *              vectors of bidiagonal matrix in VT
									// *              (Workspace: need M+BDSPAC)
									// *
									_dlaset.Run("F", N, N, ZERO, ZERO, ref VT, offset_vt
                                                     , LDVT);
									_dbdsdc.Run("L", "I", M, ref S, offset_s, ref WORK, IE + o_work, ref U, offset_u
                                                     , LDU, ref VT, offset_vt, LDVT, ref DUM, offset_dum, ref IDUM, offset_idum, ref WORK, NWORK + o_work
                                                     , ref IWORK, offset_iwork, ref INFO);
                                    // *
                                    // *              Set the right corner of VT to identity matrix
                                    // *
                                    if (N > M)
                                    {
										_dlaset.Run("F", N - M, N - M, ZERO, ONE, ref VT, M + 1+(M + 1) * LDVT + o_vt
                                                         , LDVT);
                                    }
									// *
									// *              Overwrite U by left singular vectors of A and VT
									// *              by right singular vectors of A
									// *              (Workspace: need 2*M+N, prefer 2*M+N*NB)
									// *
									_dormbr.Run("Q", "L", "N", M, M, N
                                                     , ref A, offset_a, LDA, WORK, ITAUQ + o_work, ref U, offset_u, LDU, ref WORK, NWORK + o_work
                                                     , LWORK - NWORK + 1, ref IERR);
									_dormbr.Run("P", "R", "T", N, N, M
                                                     , ref A, offset_a, LDA, WORK, ITAUP + o_work, ref VT, offset_vt, LDVT, ref WORK, NWORK + o_work
                                                     , LWORK - NWORK + 1, ref IERR);
                                }
                            }
                        }
                    }
                    // *
                }
                // *
            }
            // *
            // *     Undo scaling if necessary
            // *
            if (ISCL == 1)
            {
                if (ANRM > BIGNUM)
                {
					_dlascl.Run("G", 0, 0, BIGNUM, ANRM, MINMN
                                     , 1, ref S, offset_s, MINMN, ref IERR);
                }
                if (ANRM < SMLNUM)
                {
					_dlascl.Run("G", 0, 0, SMLNUM, ANRM, MINMN
                                     , 1, ref S, offset_s, MINMN, ref IERR);
                }
            }
            // *
            // *     Return optimal workspace in WORK(1)
            // *
            WORK[1 + o_work] = MAXWRK;
            // *
            return;
            // *
            // *     End of DGESDD
            // *

            #endregion

        }
    }
}
