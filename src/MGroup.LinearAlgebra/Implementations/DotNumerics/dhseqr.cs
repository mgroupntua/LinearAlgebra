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
    /// -- LAPACK driver routine (version 3.1) --
    /// Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
    /// November 2006
    /// Purpose
    /// =======
    /// 
    /// DHSEQR computes the eigenvalues of a Hessenberg matrix H
    /// and, optionally, the matrices T and Z from the Schur decomposition
    /// H = Z T Z**T, where T is an upper quasi-triangular matrix (the
    /// Schur form), and Z is the orthogonal matrix of Schur vectors.
    /// 
    /// Optionally Z may be postmultiplied into an input orthogonal
    /// matrix Q so that this routine can give the Schur factorization
    /// of a matrix A which has been reduced to the Hessenberg form H
    /// by the orthogonal matrix Q:  A = Q*H*Q**T = (QZ)*T*(QZ)**T.
    /// 
    ///</summary>
    public class DHSEQR
    {
    

        #region Dependencies
        
        ILAENV _ilaenv; LSAME _lsame; DLACPY _dlacpy; DLAHQR _dlahqr; DLAQR0 _dlaqr0; DLASET _dlaset; XERBLA _xerbla; 

        #endregion


        #region Variables
        
        const int NTINY = 11; const int NL = 49; const double ZERO = 0.0E0; const double ONE = 1.0E0; 
        double[] HL = new double[NL * NL];double[] WORKL = new double[NL]; 

        #endregion

        public DHSEQR(ILAENV ilaenv, LSAME lsame, DLACPY dlacpy, DLAHQR dlahqr, DLAQR0 dlaqr0, DLASET dlaset, XERBLA xerbla)
        {
    

            #region Set Dependencies
            
            _ilaenv = ilaenv; _lsame = lsame; _dlacpy = dlacpy; _dlahqr = dlahqr; _dlaqr0 = dlaqr0; 
            _dlaset = dlaset;_xerbla = xerbla; 

            #endregion

        }
    
        public DHSEQR()
        {
    

            #region Dependencies (Initialization)
            
            var ieeeck = new IEEECK();
            var iparmq = new IPARMQ();
            var lsame = new LSAME();
            var dlamc3 = new DLAMC3();
            var dcopy = new DCOPY();
            var dlabad = new DLABAD();
            var dlapy2 = new DLAPY2();
            var dnrm2 = new DNRM2();
            var dscal = new DSCAL();
            var drot = new DROT();
            var daxpy = new DAXPY();
            var xerbla = new XERBLA();
            var dlassq = new DLASSQ();
            var idamax = new IDAMAX();
            var dswap = new DSWAP();
            var dlaqr1 = new DLAQR1();
            var ilaenv = new ILAENV(ieeeck, iparmq);
            var dlacpy = new DLACPY(lsame);
            var dlamc1 = new DLAMC1(dlamc3);
            var dlamc4 = new DLAMC4(dlamc3);
            var dlamc5 = new DLAMC5(dlamc3);
            var dlamc2 = new DLAMC2(dlamc3, dlamc1, dlamc4, dlamc5);
            var dlamch = new DLAMCH(lsame, dlamc2);
            var dlanv2 = new DLANV2(dlamch, dlapy2);
            var dlarfg = new DLARFG(dlamch, dlapy2, dnrm2, dscal);
            var dlahqr = new DLAHQR(dlamch, dcopy, dlabad, dlanv2, dlarfg, drot);
            var dgemv = new DGEMV(lsame, xerbla);
            var dger = new DGER(xerbla);
            var dlarf = new DLARF(dgemv, dger, lsame);
            var dgehd2 = new DGEHD2(dlarf, dlarfg, xerbla);
            var dgemm = new DGEMM(lsame, xerbla);
            var dtrmm = new DTRMM(lsame, xerbla);
            var dtrmv = new DTRMV(lsame, xerbla);
            var dlahr2 = new DLAHR2(daxpy, dcopy, dgemm, dgemv, dlacpy, dlarfg, dscal, dtrmm, dtrmv);
            var dlarfb = new DLARFB(lsame, dcopy, dgemm, dtrmm);
            var dgehrd = new DGEHRD(daxpy, dgehd2, dgemm, dlahr2, dlarfb, dtrmm, xerbla, ilaenv);
            var dlaset = new DLASET(lsame);
            var dlarft = new DLARFT(dgemv, dtrmv, lsame);
            var dorg2r = new DORG2R(dlarf, dscal, xerbla);
            var dorgqr = new DORGQR(dlarfb, dlarft, dorg2r, xerbla, ilaenv);
            var dorghr = new DORGHR(dorgqr, xerbla, ilaenv);
            var dlange = new DLANGE(dlassq, lsame);
            var dlarfx = new DLARFX(lsame, dgemv, dger);
            var dlartg = new DLARTG(dlamch);
            var dlasy2 = new DLASY2(idamax, dlamch, dcopy, dswap);
            var dlaexc = new DLAEXC(dlamch, dlange, dlacpy, dlanv2, dlarfg, dlarfx, dlartg, dlasy2, drot);
            var dtrexc = new DTREXC(lsame, dlaexc, xerbla);
            var dlaqr2 = new DLAQR2(dlamch, dcopy, dgehrd, dgemm, dlabad, dlacpy, dlahqr, dlanv2, dlarf, dlarfg
                                       , dlaset, dorghr, dtrexc);
            var dlaqr5 = new DLAQR5(dlamch, dgemm, dlabad, dlacpy, dlaqr1, dlarfg, dlaset, dtrmm);
            var dlaqr4 = new DLAQR4(ilaenv, dlacpy, dlahqr, dlanv2, dlaqr2, dlaqr5);
            var dlaqr3 = new DLAQR3(dlamch, ilaenv, dcopy, dgehrd, dgemm, dlabad, dlacpy, dlahqr, dlanv2, dlaqr4
                                       , dlarf, dlarfg, dlaset, dorghr, dtrexc);
            var dlaqr0 = new DLAQR0(ilaenv, dlacpy, dlahqr, dlanv2, dlaqr3, dlaqr4, dlaqr5);

            #endregion


            #region Set Dependencies
            
            _ilaenv = ilaenv; _lsame = lsame; _dlacpy = dlacpy; _dlahqr = dlahqr; _dlaqr0 = dlaqr0; 
            _dlaset = dlaset;_xerbla = xerbla; 

            #endregion

        }
        /// <summary>
        /// Purpose
        /// =======
        /// 
        /// DHSEQR computes the eigenvalues of a Hessenberg matrix H
        /// and, optionally, the matrices T and Z from the Schur decomposition
        /// H = Z T Z**T, where T is an upper quasi-triangular matrix (the
        /// Schur form), and Z is the orthogonal matrix of Schur vectors.
        /// 
        /// Optionally Z may be postmultiplied into an input orthogonal
        /// matrix Q so that this routine can give the Schur factorization
        /// of a matrix A which has been reduced to the Hessenberg form H
        /// by the orthogonal matrix Q:  A = Q*H*Q**T = (QZ)*T*(QZ)**T.
        /// 
        ///</summary>
        /// <param name="JOB">
        /// (input) CHARACTER*1
        /// = 'E':  compute eigenvalues only;
        /// = 'S':  compute eigenvalues and the Schur form T.
        ///</param>
        /// <param name="COMPZ">
        /// (input) CHARACTER*1
        /// = 'N':  no Schur vectors are computed;
        /// = 'I':  Z is initialized to the unit matrix and the matrix Z
        /// of Schur vectors of H is returned;
        /// = 'V':  Z must contain an orthogonal matrix Q on entry, and
        /// the product Q*Z is returned.
        ///</param>
        /// <param name="N">
        /// (input) INTEGER
        /// The order of the matrix H.  N .GE. 0.
        ///</param>
        /// <param name="ILO">
        /// (input) INTEGER
        ///</param>
        /// <param name="IHI">
        /// (input) INTEGER
        /// It is assumed that H is already upper triangular in rows
        /// and columns 1:ILO-1 and IHI+1:N. ILO and IHI are normally
        /// set by a previous call to DGEBAL, and then passed to DGEHRD
        /// when the matrix output by DGEBAL is reduced to Hessenberg
        /// form. Otherwise ILO and IHI should be set to 1 and N
        /// respectively.  If N.GT.0, then 1.LE.ILO.LE.IHI.LE.N.
        /// If N = 0, then ILO = 1 and IHI = 0.
        ///</param>
        /// <param name="H">
        /// (input/output) DOUBLE PRECISION array, dimension (LDH,N)
        /// On entry, the upper Hessenberg matrix H.
        /// On exit, if INFO = 0 and JOB = 'S', then H contains the
        /// upper quasi-triangular matrix T from the Schur decomposition
        /// (the Schur form); 2-by-2 diagonal blocks (corresponding to
        /// complex conjugate pairs of eigenvalues) are returned in
        /// standard form, with H(i,i) = H(i+1,i+1) and
        /// H(i+1,i)*H(i,i+1).LT.0. If INFO = 0 and JOB = 'E', the
        /// contents of H are unspecified on exit.  (The output value of
        /// H when INFO.GT.0 is given under the description of INFO
        /// below.)
        /// 
        /// Unlike earlier versions of DHSEQR, this subroutine may
        /// explicitly H(i,j) = 0 for i.GT.j and j = 1, 2, ... ILO-1
        /// or j = IHI+1, IHI+2, ... N.
        ///</param>
        /// <param name="LDH">
        /// (input) INTEGER
        /// The leading dimension of the array H. LDH .GE. max(1,N).
        ///</param>
        /// <param name="WR">
        /// (output) DOUBLE PRECISION array, dimension (N)
        ///</param>
        /// <param name="WI">
        /// (output) DOUBLE PRECISION array, dimension (N)
        /// The real and imaginary parts, respectively, of the computed
        /// eigenvalues. If two eigenvalues are computed as a complex
        /// conjugate pair, they are stored in consecutive elements of
        /// WR and WI, say the i-th and (i+1)th, with WI(i) .GT. 0 and
        /// WI(i+1) .LT. 0. If JOB = 'S', the eigenvalues are stored in
        /// the same order as on the diagonal of the Schur form returned
        /// in H, with WR(i) = H(i,i) and, if H(i:i+1,i:i+1) is a 2-by-2
        /// diagonal block, WI(i) = sqrt(-H(i+1,i)*H(i,i+1)) and
        /// WI(i+1) = -WI(i).
        ///</param>
        /// <param name="Z">
        /// (input/output) DOUBLE PRECISION array, dimension (LDZ,N)
        /// If COMPZ = 'N', Z is not referenced.
        /// If COMPZ = 'I', on entry Z need not be set and on exit,
        /// if INFO = 0, Z contains the orthogonal matrix Z of the Schur
        /// vectors of H.  If COMPZ = 'V', on entry Z must contain an
        /// N-by-N matrix Q, which is assumed to be equal to the unit
        /// matrix except for the submatrix Z(ILO:IHI,ILO:IHI). On exit,
        /// if INFO = 0, Z contains Q*Z.
        /// Normally Q is the orthogonal matrix generated by DORGHR
        /// after the call to DGEHRD which formed the Hessenberg matrix
        /// H. (The output value of Z when INFO.GT.0 is given under
        /// the description of INFO below.)
        ///</param>
        /// <param name="LDZ">
        /// (input) INTEGER
        /// The leading dimension of the array Z.  if COMPZ = 'I' or
        /// COMPZ = 'V', then LDZ.GE.MAX(1,N).  Otherwize, LDZ.GE.1.
        ///</param>
        /// <param name="WORK">
        /// (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
        /// On exit, if INFO = 0, WORK(1) returns an estimate of
        /// the optimal value for LWORK.
        ///</param>
        /// <param name="LWORK">
        /// (input) INTEGER
        /// The dimension of the array WORK.  LWORK .GE. max(1,N)
        /// is sufficient, but LWORK typically as large as 6*N may
        /// be required for optimal performance.  A workspace query
        /// to determine the optimal workspace size is recommended.
        /// 
        /// If LWORK = -1, then DHSEQR does a workspace query.
        /// In this case, DHSEQR checks the input parameters and
        /// estimates the optimal workspace size for the given
        /// values of N, ILO and IHI.  The estimate is returned
        /// in WORK(1).  No error message related to LWORK is
        /// issued by XERBLA.  Neither H nor Z are accessed.
        /// 
        ///</param>
        /// <param name="INFO">
        /// (output) INTEGER
        /// =  0:  successful exit
        /// .LT. 0:  if INFO = -i, the i-th argument had an illegal
        /// value
        /// .GT. 0:  if INFO = i, DHSEQR failed to compute all of
        /// the eigenvalues.  Elements 1:ilo-1 and i+1:n of WR
        /// and WI contain those eigenvalues which have been
        /// successfully computed.  (Failures are rare.)
        /// 
        /// If INFO .GT. 0 and JOB = 'E', then on exit, the
        /// remaining unconverged eigenvalues are the eigen-
        /// values of the upper Hessenberg matrix rows and
        /// columns ILO through INFO of the final, output
        /// value of H.
        /// 
        /// If INFO .GT. 0 and JOB   = 'S', then on exit
        /// 
        /// (*)  (initial value of H)*U  = U*(final value of H)
        /// 
        /// where U is an orthogonal matrix.  The final
        /// value of H is upper Hessenberg and quasi-triangular
        /// in rows and columns INFO+1 through IHI.
        /// 
        /// If INFO .GT. 0 and COMPZ = 'V', then on exit
        /// 
        /// (final value of Z)  =  (initial value of Z)*U
        /// 
        /// where U is the orthogonal matrix in (*) (regard-
        /// less of the value of JOB.)
        /// 
        /// If INFO .GT. 0 and COMPZ = 'I', then on exit
        /// (final value of Z)  = U
        /// where U is the orthogonal matrix in (*) (regard-
        /// less of the value of JOB.)
        /// 
        /// If INFO .GT. 0 and COMPZ = 'N', then Z is not
        /// accessed.
        ///</param>
        public void Run(string JOB, string COMPZ, int N, int ILO, int IHI, ref double[] H, int offset_h
                         , int LDH, ref double[] WR, int offset_wr, ref double[] WI, int offset_wi, ref double[] Z, int offset_z, int LDZ, ref double[] WORK, int offset_work
                         , int LWORK, ref int INFO)
        {

            #region Variables
            
            var offset_hl = 0; var o_hl = -1 - NL; var offset_workl = 0; var I = 0; var KBOT = 0; var NMIN = 0; 
            var INITZ = false;var LQUERY = false; var WANTT = false; var WANTZ = false; 

            #endregion


            #region Array Index Correction
            
             var o_h = -1 - LDH + offset_h;  var o_wr = -1 + offset_wr;  var o_wi = -1 + offset_wi; 
             var o_z = -1 - LDZ + offset_z; var o_work = -1 + offset_work; 

            #endregion


            #region Strings
            
            JOB = JOB.Substring(0, 1);  COMPZ = COMPZ.Substring(0, 1);  

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
            // *     Purpose
            // *     =======
            // *
            // *     DHSEQR computes the eigenvalues of a Hessenberg matrix H
            // *     and, optionally, the matrices T and Z from the Schur decomposition
            // *     H = Z T Z**T, where T is an upper quasi-triangular matrix (the
            // *     Schur form), and Z is the orthogonal matrix of Schur vectors.
            // *
            // *     Optionally Z may be postmultiplied into an input orthogonal
            // *     matrix Q so that this routine can give the Schur factorization
            // *     of a matrix A which has been reduced to the Hessenberg form H
            // *     by the orthogonal matrix Q:  A = Q*H*Q**T = (QZ)*T*(QZ)**T.
            // *
            // *     Arguments
            // *     =========
            // *
            // *     JOB   (input) CHARACTER*1
            // *           = 'E':  compute eigenvalues only;
            // *           = 'S':  compute eigenvalues and the Schur form T.
            // *
            // *     COMPZ (input) CHARACTER*1
            // *           = 'N':  no Schur vectors are computed;
            // *           = 'I':  Z is initialized to the unit matrix and the matrix Z
            // *                   of Schur vectors of H is returned;
            // *           = 'V':  Z must contain an orthogonal matrix Q on entry, and
            // *                   the product Q*Z is returned.
            // *
            // *     N     (input) INTEGER
            // *           The order of the matrix H.  N .GE. 0.
            // *
            // *     ILO   (input) INTEGER
            // *     IHI   (input) INTEGER
            // *           It is assumed that H is already upper triangular in rows
            // *           and columns 1:ILO-1 and IHI+1:N. ILO and IHI are normally
            // *           set by a previous call to DGEBAL, and then passed to DGEHRD
            // *           when the matrix output by DGEBAL is reduced to Hessenberg
            // *           form. Otherwise ILO and IHI should be set to 1 and N
            // *           respectively.  If N.GT.0, then 1.LE.ILO.LE.IHI.LE.N.
            // *           If N = 0, then ILO = 1 and IHI = 0.
            // *
            // *     H     (input/output) DOUBLE PRECISION array, dimension (LDH,N)
            // *           On entry, the upper Hessenberg matrix H.
            // *           On exit, if INFO = 0 and JOB = 'S', then H contains the
            // *           upper quasi-triangular matrix T from the Schur decomposition
            // *           (the Schur form); 2-by-2 diagonal blocks (corresponding to
            // *           complex conjugate pairs of eigenvalues) are returned in
            // *           standard form, with H(i,i) = H(i+1,i+1) and
            // *           H(i+1,i)*H(i,i+1).LT.0. If INFO = 0 and JOB = 'E', the
            // *           contents of H are unspecified on exit.  (The output value of
            // *           H when INFO.GT.0 is given under the description of INFO
            // *           below.)
            // *
            // *           Unlike earlier versions of DHSEQR, this subroutine may
            // *           explicitly H(i,j) = 0 for i.GT.j and j = 1, 2, ... ILO-1
            // *           or j = IHI+1, IHI+2, ... N.
            // *
            // *     LDH   (input) INTEGER
            // *           The leading dimension of the array H. LDH .GE. max(1,N).
            // *
            // *     WR    (output) DOUBLE PRECISION array, dimension (N)
            // *     WI    (output) DOUBLE PRECISION array, dimension (N)
            // *           The real and imaginary parts, respectively, of the computed
            // *           eigenvalues. If two eigenvalues are computed as a complex
            // *           conjugate pair, they are stored in consecutive elements of
            // *           WR and WI, say the i-th and (i+1)th, with WI(i) .GT. 0 and
            // *           WI(i+1) .LT. 0. If JOB = 'S', the eigenvalues are stored in
            // *           the same order as on the diagonal of the Schur form returned
            // *           in H, with WR(i) = H(i,i) and, if H(i:i+1,i:i+1) is a 2-by-2
            // *           diagonal block, WI(i) = sqrt(-H(i+1,i)*H(i,i+1)) and
            // *           WI(i+1) = -WI(i).
            // *
            // *     Z     (input/output) DOUBLE PRECISION array, dimension (LDZ,N)
            // *           If COMPZ = 'N', Z is not referenced.
            // *           If COMPZ = 'I', on entry Z need not be set and on exit,
            // *           if INFO = 0, Z contains the orthogonal matrix Z of the Schur
            // *           vectors of H.  If COMPZ = 'V', on entry Z must contain an
            // *           N-by-N matrix Q, which is assumed to be equal to the unit
            // *           matrix except for the submatrix Z(ILO:IHI,ILO:IHI). On exit,
            // *           if INFO = 0, Z contains Q*Z.
            // *           Normally Q is the orthogonal matrix generated by DORGHR
            // *           after the call to DGEHRD which formed the Hessenberg matrix
            // *           H. (The output value of Z when INFO.GT.0 is given under
            // *           the description of INFO below.)
            // *
            // *     LDZ   (input) INTEGER
            // *           The leading dimension of the array Z.  if COMPZ = 'I' or
            // *           COMPZ = 'V', then LDZ.GE.MAX(1,N).  Otherwize, LDZ.GE.1.
            // *
            // *     WORK  (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
            // *           On exit, if INFO = 0, WORK(1) returns an estimate of
            // *           the optimal value for LWORK.
            // *
            // *     LWORK (input) INTEGER
            // *           The dimension of the array WORK.  LWORK .GE. max(1,N)
            // *           is sufficient, but LWORK typically as large as 6*N may
            // *           be required for optimal performance.  A workspace query
            // *           to determine the optimal workspace size is recommended.
            // *
            // *           If LWORK = -1, then DHSEQR does a workspace query.
            // *           In this case, DHSEQR checks the input parameters and
            // *           estimates the optimal workspace size for the given
            // *           values of N, ILO and IHI.  The estimate is returned
            // *           in WORK(1).  No error message related to LWORK is
            // *           issued by XERBLA.  Neither H nor Z are accessed.
            // *
            // *
            // *     INFO  (output) INTEGER
            // *             =  0:  successful exit
            // *           .LT. 0:  if INFO = -i, the i-th argument had an illegal
            // *                    value
            // *           .GT. 0:  if INFO = i, DHSEQR failed to compute all of
            // *                the eigenvalues.  Elements 1:ilo-1 and i+1:n of WR
            // *                and WI contain those eigenvalues which have been
            // *                successfully computed.  (Failures are rare.)
            // *
            // *                If INFO .GT. 0 and JOB = 'E', then on exit, the
            // *                remaining unconverged eigenvalues are the eigen-
            // *                values of the upper Hessenberg matrix rows and
            // *                columns ILO through INFO of the final, output
            // *                value of H.
            // *
            // *                If INFO .GT. 0 and JOB   = 'S', then on exit
            // *
            // *           (*)  (initial value of H)*U  = U*(final value of H)
            // *
            // *                where U is an orthogonal matrix.  The final
            // *                value of H is upper Hessenberg and quasi-triangular
            // *                in rows and columns INFO+1 through IHI.
            // *
            // *                If INFO .GT. 0 and COMPZ = 'V', then on exit
            // *
            // *                  (final value of Z)  =  (initial value of Z)*U
            // *
            // *                where U is the orthogonal matrix in (*) (regard-
            // *                less of the value of JOB.)
            // *
            // *                If INFO .GT. 0 and COMPZ = 'I', then on exit
            // *                      (final value of Z)  = U
            // *                where U is the orthogonal matrix in (*) (regard-
            // *                less of the value of JOB.)
            // *
            // *                If INFO .GT. 0 and COMPZ = 'N', then Z is not
            // *                accessed.
            // *
            // *     ================================================================
            // *             Default values supplied by
            // *             ILAENV(ISPEC,'DHSEQR',JOB(:1)//COMPZ(:1),N,ILO,IHI,LWORK).
            // *             It is suggested that these defaults be adjusted in order
            // *             to attain best performance in each particular
            // *             computational environment.
            // *
            // *            ISPEC=1:  The DLAHQR vs DLAQR0 crossover point.
            // *                      Default: 75. (Must be at least 11.)
            // *
            // *            ISPEC=2:  Recommended deflation window size.
            // *                      This depends on ILO, IHI and NS.  NS is the
            // *                      number of simultaneous shifts returned
            // *                      by ILAENV(ISPEC=4).  (See ISPEC=4 below.)
            // *                      The default for (IHI-ILO+1).LE.500 is NS.
            // *                      The default for (IHI-ILO+1).GT.500 is 3*NS/2.
            // *
            // *            ISPEC=3:  Nibble crossover point. (See ILAENV for
            // *                      details.)  Default: 14% of deflation window
            // *                      size.
            // *
            // *            ISPEC=4:  Number of simultaneous shifts, NS, in
            // *                      a multi-shift QR iteration.
            // *
            // *                      If IHI-ILO+1 is ...
            // *
            // *                      greater than      ...but less    ... the
            // *                      or equal to ...      than        default is
            // *
            // *                           1               30          NS -   2(+)
            // *                          30               60          NS -   4(+)
            // *                          60              150          NS =  10(+)
            // *                         150              590          NS =  **
            // *                         590             3000          NS =  64
            // *                        3000             6000          NS = 128
            // *                        6000             infinity      NS = 256
            // *
            // *                  (+)  By default some or all matrices of this order 
            // *                       are passed to the implicit double shift routine
            // *                       DLAHQR and NS is ignored.  See ISPEC=1 above 
            // *                       and comments in IPARM for details.
            // *
            // *                       The asterisks (**) indicate an ad-hoc
            // *                       function of N increasing from 10 to 64.
            // *
            // *            ISPEC=5:  Select structured matrix multiply.
            // *                      (See ILAENV for details.) Default: 3.
            // *
            // *     ================================================================
            // *     Based on contributions by
            // *        Karen Braman and Ralph Byers, Department of Mathematics,
            // *        University of Kansas, USA
            // *
            // *     ================================================================
            // *     References:
            // *       K. Braman, R. Byers and R. Mathias, The Multi-Shift QR
            // *       Algorithm Part I: Maintaining Well Focused Shifts, and Level 3
            // *       Performance, SIAM Journal of Matrix Analysis, volume 23, pages
            // *       929--947, 2002.
            // *
            // *       K. Braman, R. Byers and R. Mathias, The Multi-Shift QR
            // *       Algorithm Part II: Aggressive Early Deflation, SIAM Journal
            // *       of Matrix Analysis, volume 23, pages 948--973, 2002.
            // *
            // *     ================================================================
            // *     .. Parameters ..
            // *
            // *     ==== Matrices of order NTINY or smaller must be processed by
            // *     .    DLAHQR because of insufficient subdiagonal scratch space.
            // *     .    (This is a hard limit.) ====
            // *
            // *     ==== NL allocates some local workspace to help small matrices
            // *     .    through a rare DLAHQR failure.  NL .GT. NTINY = 11 is
            // *     .    required and NL .LE. NMIN = ILAENV(ISPEC=1,...) is recom-
            // *     .    mended.  (The default value of NMIN is 75.)  Using NL = 49
            // *     .    allows up to six simultaneous shifts and a 16-by-16
            // *     .    deflation window.  ====
            // *
            // *     ..
            // *     .. Local Arrays ..
            // *     ..
            // *     .. Local Scalars ..
            // *     ..
            // *     .. External Functions ..
            // *     ..
            // *     .. External Subroutines ..
            // *     ..
            // *     .. Intrinsic Functions ..
            //      INTRINSIC          DBLE, MAX, MIN;
            // *     ..
            // *     .. Executable Statements ..
            // *
            // *     ==== Decode and check the input parameters. ====
            // *

            #endregion


            #region Body
            
            WANTT = _lsame.Run(JOB, "S");
            INITZ = _lsame.Run(COMPZ, "I");
            WANTZ = INITZ || _lsame.Run(COMPZ, "V");
            WORK[1 + o_work] = Convert.ToDouble(Math.Max(1, N));
            LQUERY = LWORK ==  - 1;
            // *
            INFO = 0;
            if (!_lsame.Run(JOB, "E") && !WANTT)
            {
                INFO =  - 1;
            }
            else
            {
                if (!_lsame.Run(COMPZ, "N") && !WANTZ)
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
                        if (ILO < 1 || ILO > Math.Max(1, N))
                        {
                            INFO =  - 4;
                        }
                        else
                        {
                            if (IHI < Math.Min(ILO, N) || IHI > N)
                            {
                                INFO =  - 5;
                            }
                            else
                            {
                                if (LDH < Math.Max(1, N))
                                {
                                    INFO =  - 7;
                                }
                                else
                                {
                                    if (LDZ < 1 || WANTZ && LDZ < Math.Max(1, N))
                                    {
                                        INFO =  - 11;
                                    }
                                    else
                                    {
                                        if (LWORK < Math.Max(1, N) && !LQUERY)
                                        {
                                            INFO =  - 13;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            // *
            if (INFO != 0)
            {
                // *
                // *        ==== Quick return in case of invalid argument. ====
                // *
                _xerbla.Run("DHSEQR",  - INFO);
                return;
                // *
            }
            else
            {
                if (N == 0)
                {
                    // *
                    // *        ==== Quick return in case N = 0; nothing to do. ====
                    // *
                    return;
                    // *
                }
                else
                {
                    if (LQUERY)
                    {
                        // *
                        // *        ==== Quick return in case of a workspace query ====
                        // *
                        _dlaqr0.Run(WANTT, WANTZ, N, ILO, IHI, ref H, offset_h
                                         , LDH, ref WR, offset_wr, ref WI, offset_wi, ILO, IHI, ref Z, offset_z
                                         , LDZ, ref WORK, offset_work, LWORK, ref INFO);
                        // *        ==== Ensure reported workspace size is backward-compatible with
                        // *        .    previous LAPACK versions. ====
                        WORK[1 + o_work] = Math.Max(Convert.ToDouble(Math.Max(1, N)), WORK[1 + o_work]);
                        return;
                        // *
                    }
                    else
                    {
                        // *
                        // *        ==== copy eigenvalues isolated by DGEBAL ====
                        // *
                        for (I = 1; I <= ILO - 1; I++)
                        {
                            WR[I + o_wr] = H[I+I * LDH + o_h];
                            WI[I + o_wi] = ZERO;
                        }
                        for (I = IHI + 1; I <= N; I++)
                        {
                            WR[I + o_wr] = H[I+I * LDH + o_h];
                            WI[I + o_wi] = ZERO;
                        }
                        // *
                        // *        ==== Initialize Z, if requested ====
                        // *
                        if (INITZ)
                        {
                            _dlaset.Run("A", N, N, ZERO, ONE, ref Z, offset_z
                                             , LDZ);
                        }
                        // *
                        // *        ==== Quick return if possible ====
                        // *
                        if (ILO == IHI)
                        {
                            WR[ILO + o_wr] = H[ILO+ILO * LDH + o_h];
                            WI[ILO + o_wi] = ZERO;
                            return;
                        }
                        // *
                        // *        ==== DLAHQR/DLAQR0 crossover point ====
                        // *
                        NMIN = _ilaenv.Run(12, "DHSEQR", FortranLib.Substring(JOB, 1, 1) + FortranLib.Substring(COMPZ, 1, 1), N, ILO, IHI, LWORK);
                        NMIN = Math.Max(NTINY, NMIN);
                        // *
                        // *        ==== DLAQR0 for big matrices; DLAHQR for small ones ====
                        // *
                        if (N > NMIN)
                        {
                            _dlaqr0.Run(WANTT, WANTZ, N, ILO, IHI, ref H, offset_h
                                             , LDH, ref WR, offset_wr, ref WI, offset_wi, ILO, IHI, ref Z, offset_z
                                             , LDZ, ref WORK, offset_work, LWORK, ref INFO);
                        }
                        else
                        {
                            // *
                            // *           ==== Small matrix ====
                            // *
                            _dlahqr.Run(WANTT, WANTZ, N, ILO, IHI, ref H, offset_h
                                             , LDH, ref WR, offset_wr, ref WI, offset_wi, ILO, IHI, ref Z, offset_z
                                             , LDZ, ref INFO);
                            // *
                            if (INFO > 0)
                            {
                                // *
                                // *              ==== A rare DLAHQR failure!  DLAQR0 sometimes succeeds
                                // *              .    when DLAHQR fails. ====
                                // *
                                KBOT = INFO;
                                // *
                                if (N >= NL)
                                {
                                    // *
                                    // *                 ==== Larger matrices have enough subdiagonal scratch
                                    // *                 .    space to call DLAQR0 directly. ====
                                    // *
                                    _dlaqr0.Run(WANTT, WANTZ, N, ILO, KBOT, ref H, offset_h
                                                     , LDH, ref WR, offset_wr, ref WI, offset_wi, ILO, IHI, ref Z, offset_z
                                                     , LDZ, ref WORK, offset_work, LWORK, ref INFO);
                                    // *
                                }
                                else
                                {
                                    // *
                                    // *                 ==== Tiny matrices don't have enough subdiagonal
                                    // *                 .    scratch space to benefit from DLAQR0.  Hence,
                                    // *                 .    tiny matrices must be copied into a larger
                                    // *                 .    array before calling DLAQR0. ====
                                    // *
                                    _dlacpy.Run("A", N, N, H, offset_h, LDH, ref HL, offset_hl
                                                     , NL);
                                    HL[N + 1+N * NL + o_hl] = ZERO;
                                    _dlaset.Run("A", NL, NL - N, ZERO, ZERO, ref HL, 1+(N + 1) * NL + o_hl
                                                     , NL);
                                    _dlaqr0.Run(WANTT, WANTZ, NL, ILO, KBOT, ref HL, offset_hl
                                                     , NL, ref WR, offset_wr, ref WI, offset_wi, ILO, IHI, ref Z, offset_z
                                                     , LDZ, ref WORKL, offset_workl, NL, ref INFO);
                                    if (WANTT || INFO != 0)
                                    {
                                        _dlacpy.Run("A", N, N, HL, offset_hl, NL, ref H, offset_h
                                                         , LDH);
                                    }
                                }
                            }
                        }
                        // *
                        // *        ==== Clear out the trash, if necessary. ====
                        // *
                        if ((WANTT || INFO != 0) && N > 2)
                        {
                            _dlaset.Run("L", N - 2, N - 2, ZERO, ZERO, ref H, 3+1 * LDH + o_h
                                             , LDH);
                        }
                        // *
                        // *        ==== Ensure reported workspace size is backward-compatible with
                        // *        .    previous LAPACK versions. ====
                        // *
                        WORK[1 + o_work] = Math.Max(Convert.ToDouble(Math.Max(1, N)), WORK[1 + o_work]);
                    }
                }
            }
            // *
            // *     ==== End of DHSEQR ====
            // *

            #endregion

        }
    }
}
