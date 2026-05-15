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
    ///</summary>
    public class DLAQR2
    {
    

        #region Dependencies
        
        DLAMCH _dlamch; DCOPY _dcopy; DGEHRD _dgehrd; DGEMM _dgemm; DLABAD _dlabad; DLACPY _dlacpy; DLAHQR _dlahqr; 
        DLANV2 _dlanv2;DLARF _dlarf; DLARFG _dlarfg; DLASET _dlaset; DORGHR _dorghr; DTREXC _dtrexc; 

        #endregion


        #region Variables
        
        const double ZERO = 0.0E0; const double ONE = 1.0E0; 

        #endregion

        public DLAQR2(DLAMCH dlamch, DCOPY dcopy, DGEHRD dgehrd, DGEMM dgemm, DLABAD dlabad, DLACPY dlacpy, DLAHQR dlahqr, DLANV2 dlanv2, DLARF dlarf, DLARFG dlarfg
                      , DLASET dlaset, DORGHR dorghr, DTREXC dtrexc)
        {
    

            #region Set Dependencies
            
            _dlamch = dlamch; _dcopy = dcopy; _dgehrd = dgehrd; _dgemm = dgemm; _dlabad = dlabad; 
            _dlacpy = dlacpy;_dlahqr = dlahqr; _dlanv2 = dlanv2; _dlarf = dlarf; _dlarfg = dlarfg; 
            _dlaset = dlaset;_dorghr = dorghr; _dtrexc = dtrexc; 

            #endregion

        }
    
        public DLAQR2()
        {
    

            #region Dependencies (Initialization)
            
            var lsame = new LSAME();
            var dlamc3 = new DLAMC3();
            var dcopy = new DCOPY();
            var daxpy = new DAXPY();
            var xerbla = new XERBLA();
            var dlapy2 = new DLAPY2();
            var dnrm2 = new DNRM2();
            var dscal = new DSCAL();
            var ieeeck = new IEEECK();
            var iparmq = new IPARMQ();
            var dlabad = new DLABAD();
            var drot = new DROT();
            var dlassq = new DLASSQ();
            var idamax = new IDAMAX();
            var dswap = new DSWAP();
            var dlamc1 = new DLAMC1(dlamc3);
            var dlamc4 = new DLAMC4(dlamc3);
            var dlamc5 = new DLAMC5(dlamc3);
            var dlamc2 = new DLAMC2(dlamc3, dlamc1, dlamc4, dlamc5);
            var dlamch = new DLAMCH(lsame, dlamc2);
            var dgemv = new DGEMV(lsame, xerbla);
            var dger = new DGER(xerbla);
            var dlarf = new DLARF(dgemv, dger, lsame);
            var dlarfg = new DLARFG(dlamch, dlapy2, dnrm2, dscal);
            var dgehd2 = new DGEHD2(dlarf, dlarfg, xerbla);
            var dgemm = new DGEMM(lsame, xerbla);
            var dlacpy = new DLACPY(lsame);
            var dtrmm = new DTRMM(lsame, xerbla);
            var dtrmv = new DTRMV(lsame, xerbla);
            var dlahr2 = new DLAHR2(daxpy, dcopy, dgemm, dgemv, dlacpy, dlarfg, dscal, dtrmm, dtrmv);
            var dlarfb = new DLARFB(lsame, dcopy, dgemm, dtrmm);
            var ilaenv = new ILAENV(ieeeck, iparmq);
            var dgehrd = new DGEHRD(daxpy, dgehd2, dgemm, dlahr2, dlarfb, dtrmm, xerbla, ilaenv);
            var dlanv2 = new DLANV2(dlamch, dlapy2);
            var dlahqr = new DLAHQR(dlamch, dcopy, dlabad, dlanv2, dlarfg, drot);
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

            #endregion


            #region Set Dependencies
            
            _dlamch = dlamch; _dcopy = dcopy; _dgehrd = dgehrd; _dgemm = dgemm; _dlabad = dlabad; 
            _dlacpy = dlacpy;_dlahqr = dlahqr; _dlanv2 = dlanv2; _dlarf = dlarf; _dlarfg = dlarfg; 
            _dlaset = dlaset;_dorghr = dorghr; _dtrexc = dtrexc; 

            #endregion

        }
        /// <param name="WANTT">
        /// (input) LOGICAL
        /// If .TRUE., then the Hessenberg matrix H is fully updated
        /// so that the quasi-triangular Schur factor may be
        /// computed (in cooperation with the calling subroutine).
        /// If .FALSE., then only enough of H is updated to preserve
        /// the eigenvalues.
        ///</param>
        /// <param name="WANTZ">
        /// (input) LOGICAL
        /// If .TRUE., then the orthogonal matrix Z is updated so
        /// so that the orthogonal Schur factor may be computed
        /// (in cooperation with the calling subroutine).
        /// If .FALSE., then Z is not referenced.
        ///</param>
        /// <param name="N">
        /// (input) INTEGER
        /// The order of the matrix H and (if WANTZ is .TRUE.) the
        /// order of the orthogonal matrix Z.
        ///</param>
        /// <param name="KTOP">
        /// (input) INTEGER
        /// It is assumed that either KTOP = 1 or H(KTOP,KTOP-1)=0.
        /// KBOT and KTOP together determine an isolated block
        /// along the diagonal of the Hessenberg matrix.
        ///</param>
        /// <param name="KBOT">
        /// (input) INTEGER
        /// It is assumed without a check that either
        /// KBOT = N or H(KBOT+1,KBOT)=0.  KBOT and KTOP together
        /// determine an isolated block along the diagonal of the
        /// Hessenberg matrix.
        ///</param>
        /// <param name="NW">
        /// (input) INTEGER
        /// Deflation window size.  1 .LE. NW .LE. (KBOT-KTOP+1).
        ///</param>
        /// <param name="H">
        /// (input/output) DOUBLE PRECISION array, dimension (LDH,N)
        /// On input the initial N-by-N section of H stores the
        /// Hessenberg matrix undergoing aggressive early deflation.
        /// On output H has been transformed by an orthogonal
        /// similarity transformation, perturbed, and the returned
        /// to Hessenberg form that (it is to be hoped) has some
        /// zero subdiagonal entries.
        ///</param>
        /// <param name="LDH">
        /// (input) integer
        /// Leading dimension of H just as declared in the calling
        /// subroutine.  N .LE. LDH
        ///</param>
        /// <param name="ILOZ">
        /// (input) INTEGER
        ///</param>
        /// <param name="IHIZ">
        /// (input) INTEGER
        /// Specify the rows of Z to which transformations must be
        /// applied if WANTZ is .TRUE.. 1 .LE. ILOZ .LE. IHIZ .LE. N.
        ///</param>
        /// <param name="Z">
        /// (input/output) DOUBLE PRECISION array, dimension (LDZ,IHI)
        /// IF WANTZ is .TRUE., then on output, the orthogonal
        /// similarity transformation mentioned above has been
        /// accumulated into Z(ILOZ:IHIZ,ILO:IHI) from the right.
        /// If WANTZ is .FALSE., then Z is unreferenced.
        ///</param>
        /// <param name="LDZ">
        /// (input) integer
        /// The leading dimension of Z just as declared in the
        /// calling subroutine.  1 .LE. LDZ.
        ///</param>
        /// <param name="NS">
        /// (output) integer
        /// The number of unconverged (ie approximate) eigenvalues
        /// returned in SR and SI that may be used as shifts by the
        /// calling subroutine.
        ///</param>
        /// <param name="ND">
        /// (output) integer
        /// The number of converged eigenvalues uncovered by this
        /// subroutine.
        ///</param>
        /// <param name="SR">
        /// (output) DOUBLE PRECISION array, dimension KBOT
        ///</param>
        /// <param name="SI">
        /// (output) DOUBLE PRECISION array, dimension KBOT
        /// On output, the real and imaginary parts of approximate
        /// eigenvalues that may be used for shifts are stored in
        /// SR(KBOT-ND-NS+1) through SR(KBOT-ND) and
        /// SI(KBOT-ND-NS+1) through SI(KBOT-ND), respectively.
        /// The real and imaginary parts of converged eigenvalues
        /// are stored in SR(KBOT-ND+1) through SR(KBOT) and
        /// SI(KBOT-ND+1) through SI(KBOT), respectively.
        ///</param>
        /// <param name="V">
        /// (workspace) DOUBLE PRECISION array, dimension (LDV,NW)
        /// An NW-by-NW work array.
        ///</param>
        /// <param name="LDV">
        /// (input) integer scalar
        /// The leading dimension of V just as declared in the
        /// calling subroutine.  NW .LE. LDV
        ///</param>
        /// <param name="NH">
        /// (input) integer scalar
        /// The number of columns of T.  NH.GE.NW.
        ///</param>
        /// <param name="T">
        /// (workspace) DOUBLE PRECISION array, dimension (LDT,NW)
        ///</param>
        /// <param name="LDT">
        /// (input) integer
        /// The leading dimension of T just as declared in the
        /// calling subroutine.  NW .LE. LDT
        ///</param>
        /// <param name="NV">
        /// (input) integer
        /// The number of rows of work array WV available for
        /// workspace.  NV.GE.NW.
        ///</param>
        /// <param name="WV">
        /// (workspace) DOUBLE PRECISION array, dimension (LDWV,NW)
        ///</param>
        /// <param name="LDWV">
        /// (input) integer
        /// The leading dimension of W just as declared in the
        /// calling subroutine.  NW .LE. LDV
        ///</param>
        /// <param name="WORK">
        /// (workspace) DOUBLE PRECISION array, dimension LWORK.
        /// On exit, WORK(1) is set to an estimate of the optimal value
        /// of LWORK for the given values of N, NW, KTOP and KBOT.
        ///</param>
        /// <param name="LWORK">
        /// (input) integer
        /// The dimension of the work array WORK.  LWORK = 2*NW
        /// suffices, but greater efficiency may result from larger
        /// values of LWORK.
        /// 
        /// If LWORK = -1, then a workspace query is assumed; DLAQR2
        /// only estimates the optimal workspace size for the given
        /// values of N, NW, KTOP and KBOT.  The estimate is returned
        /// in WORK(1).  No error message related to LWORK is issued
        /// by XERBLA.  Neither H nor Z are accessed.
        ///</param>
        public void Run(bool WANTT, bool WANTZ, int N, int KTOP, int KBOT, int NW
                         , ref double[] H, int offset_h, int LDH, int ILOZ, int IHIZ, ref double[] Z, int offset_z, int LDZ
                         , ref int NS, ref int ND, ref double[] SR, int offset_sr, ref double[] SI, int offset_si, ref double[] V, int offset_v, int LDV
                         , int NH, ref double[] T, int offset_t, int LDT, int NV, ref double[] WV, int offset_wv, int LDWV
                         , ref double[] WORK, int offset_work, int LWORK)
        {

            #region Variables
            
            double AA = 0; double BB = 0; double BETA = 0; double CC = 0; double CS = 0; double DD = 0; double EVI = 0; 
            double EVK = 0;double FOO = 0; double S = 0; double SAFMAX = 0; double SAFMIN = 0; double SMLNUM = 0; double SN = 0; 
            double TAU = 0;double ULP = 0; var I = 0; var IFST = 0; var ILST = 0; var INFO = 0; var INFQR = 0; var J = 0; 
            var JW = 0;var K = 0; var KCOL = 0; var KEND = 0; var KLN = 0; var KROW = 0; var KWTOP = 0; var LTOP = 0; 
            var LWK1 = 0;var LWK2 = 0; var LWKOPT = 0; var BULGE = false; var SORTED = false; 

            #endregion


            #region Array Index Correction
            
             var o_h = -1 - LDH + offset_h;  var o_z = -1 - LDZ + offset_z;  var o_sr = -1 + offset_sr; 
             var o_si = -1 + offset_si; var o_v = -1 - LDV + offset_v;  var o_t = -1 - LDT + offset_t; 
             var o_wv = -1 - LDWV + offset_wv; var o_work = -1 + offset_work; 

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
            // *     This subroutine is identical to DLAQR3 except that it avoids
            // *     recursion by calling DLAHQR instead of DLAQR4.
            // *
            // *
            // *     ******************************************************************
            // *     Aggressive early deflation:
            // *
            // *     This subroutine accepts as input an upper Hessenberg matrix
            // *     H and performs an orthogonal similarity transformation
            // *     designed to detect and deflate fully converged eigenvalues from
            // *     a trailing principal submatrix.  On output H has been over-
            // *     written by a new Hessenberg matrix that is a perturbation of
            // *     an orthogonal similarity transformation of H.  It is to be
            // *     hoped that the final version of H has many zero subdiagonal
            // *     entries.
            // *
            // *     ******************************************************************
            // *     WANTT   (input) LOGICAL
            // *          If .TRUE., then the Hessenberg matrix H is fully updated
            // *          so that the quasi-triangular Schur factor may be
            // *          computed (in cooperation with the calling subroutine).
            // *          If .FALSE., then only enough of H is updated to preserve
            // *          the eigenvalues.
            // *
            // *     WANTZ   (input) LOGICAL
            // *          If .TRUE., then the orthogonal matrix Z is updated so
            // *          so that the orthogonal Schur factor may be computed
            // *          (in cooperation with the calling subroutine).
            // *          If .FALSE., then Z is not referenced.
            // *
            // *     N       (input) INTEGER
            // *          The order of the matrix H and (if WANTZ is .TRUE.) the
            // *          order of the orthogonal matrix Z.
            // *
            // *     KTOP    (input) INTEGER
            // *          It is assumed that either KTOP = 1 or H(KTOP,KTOP-1)=0.
            // *          KBOT and KTOP together determine an isolated block
            // *          along the diagonal of the Hessenberg matrix.
            // *
            // *     KBOT    (input) INTEGER
            // *          It is assumed without a check that either
            // *          KBOT = N or H(KBOT+1,KBOT)=0.  KBOT and KTOP together
            // *          determine an isolated block along the diagonal of the
            // *          Hessenberg matrix.
            // *
            // *     NW      (input) INTEGER
            // *          Deflation window size.  1 .LE. NW .LE. (KBOT-KTOP+1).
            // *
            // *     H       (input/output) DOUBLE PRECISION array, dimension (LDH,N)
            // *          On input the initial N-by-N section of H stores the
            // *          Hessenberg matrix undergoing aggressive early deflation.
            // *          On output H has been transformed by an orthogonal
            // *          similarity transformation, perturbed, and the returned
            // *          to Hessenberg form that (it is to be hoped) has some
            // *          zero subdiagonal entries.
            // *
            // *     LDH     (input) integer
            // *          Leading dimension of H just as declared in the calling
            // *          subroutine.  N .LE. LDH
            // *
            // *     ILOZ    (input) INTEGER
            // *     IHIZ    (input) INTEGER
            // *          Specify the rows of Z to which transformations must be
            // *          applied if WANTZ is .TRUE.. 1 .LE. ILOZ .LE. IHIZ .LE. N.
            // *
            // *     Z       (input/output) DOUBLE PRECISION array, dimension (LDZ,IHI)
            // *          IF WANTZ is .TRUE., then on output, the orthogonal
            // *          similarity transformation mentioned above has been
            // *          accumulated into Z(ILOZ:IHIZ,ILO:IHI) from the right.
            // *          If WANTZ is .FALSE., then Z is unreferenced.
            // *
            // *     LDZ     (input) integer
            // *          The leading dimension of Z just as declared in the
            // *          calling subroutine.  1 .LE. LDZ.
            // *
            // *     NS      (output) integer
            // *          The number of unconverged (ie approximate) eigenvalues
            // *          returned in SR and SI that may be used as shifts by the
            // *          calling subroutine.
            // *
            // *     ND      (output) integer
            // *          The number of converged eigenvalues uncovered by this
            // *          subroutine.
            // *
            // *     SR      (output) DOUBLE PRECISION array, dimension KBOT
            // *     SI      (output) DOUBLE PRECISION array, dimension KBOT
            // *          On output, the real and imaginary parts of approximate
            // *          eigenvalues that may be used for shifts are stored in
            // *          SR(KBOT-ND-NS+1) through SR(KBOT-ND) and
            // *          SI(KBOT-ND-NS+1) through SI(KBOT-ND), respectively.
            // *          The real and imaginary parts of converged eigenvalues
            // *          are stored in SR(KBOT-ND+1) through SR(KBOT) and
            // *          SI(KBOT-ND+1) through SI(KBOT), respectively.
            // *
            // *     V       (workspace) DOUBLE PRECISION array, dimension (LDV,NW)
            // *          An NW-by-NW work array.
            // *
            // *     LDV     (input) integer scalar
            // *          The leading dimension of V just as declared in the
            // *          calling subroutine.  NW .LE. LDV
            // *
            // *     NH      (input) integer scalar
            // *          The number of columns of T.  NH.GE.NW.
            // *
            // *     T       (workspace) DOUBLE PRECISION array, dimension (LDT,NW)
            // *
            // *     LDT     (input) integer
            // *          The leading dimension of T just as declared in the
            // *          calling subroutine.  NW .LE. LDT
            // *
            // *     NV      (input) integer
            // *          The number of rows of work array WV available for
            // *          workspace.  NV.GE.NW.
            // *
            // *     WV      (workspace) DOUBLE PRECISION array, dimension (LDWV,NW)
            // *
            // *     LDWV    (input) integer
            // *          The leading dimension of W just as declared in the
            // *          calling subroutine.  NW .LE. LDV
            // *
            // *     WORK    (workspace) DOUBLE PRECISION array, dimension LWORK.
            // *          On exit, WORK(1) is set to an estimate of the optimal value
            // *          of LWORK for the given values of N, NW, KTOP and KBOT.
            // *
            // *     LWORK   (input) integer
            // *          The dimension of the work array WORK.  LWORK = 2*NW
            // *          suffices, but greater efficiency may result from larger
            // *          values of LWORK.
            // *
            // *          If LWORK = -1, then a workspace query is assumed; DLAQR2
            // *          only estimates the optimal workspace size for the given
            // *          values of N, NW, KTOP and KBOT.  The estimate is returned
            // *          in WORK(1).  No error message related to LWORK is issued
            // *          by XERBLA.  Neither H nor Z are accessed.
            // *
            // *     ================================================================
            // *     Based on contributions by
            // *        Karen Braman and Ralph Byers, Department of Mathematics,
            // *        University of Kansas, USA
            // *
            // *     ================================================================
            // *     .. Parameters ..
            // *     ..
            // *     .. Local Scalars ..
            // *     ..
            // *     .. External Functions ..
            // *     ..
            // *     .. External Subroutines ..
            // *     ..
            // *     .. Intrinsic Functions ..
            //      INTRINSIC          ABS, DBLE, INT, MAX, MIN, SQRT;
            // *     ..
            // *     .. Executable Statements ..
            // *
            // *     ==== Estimate optimal workspace. ====
            // *

            #endregion


            #region Body
            
            JW = Math.Min(NW, KBOT - KTOP + 1);
            if (JW <= 2)
            {
                LWKOPT = 1;
            }
            else
            {
                // *
                // *        ==== Workspace query call to DGEHRD ====
                // *
                _dgehrd.Run(JW, 1, JW - 1, ref T, offset_t, LDT, ref WORK, offset_work
                                 , ref WORK, offset_work,  - 1, ref INFO);
                LWK1 = Convert.ToInt32(Math.Truncate(WORK[1 + o_work]));
                // *
                // *        ==== Workspace query call to DORGHR ====
                // *
                _dorghr.Run(JW, 1, JW - 1, ref T, offset_t, LDT, WORK, offset_work
                                 , ref WORK, offset_work,  - 1, ref INFO);
                LWK2 = Convert.ToInt32(Math.Truncate(WORK[1 + o_work]));
                // *
                // *        ==== Optimal workspace ====
                // *
                LWKOPT = JW + Math.Max(LWK1, LWK2);
            }
            // *
            // *     ==== Quick return in case of workspace query. ====
            // *
            if (LWORK ==  - 1)
            {
                WORK[1 + o_work] = Convert.ToDouble(LWKOPT);
                return;
            }
            // *
            // *     ==== Nothing to do ...
            // *     ... for an empty active block ... ====
            NS = 0;
            ND = 0;
            if (KTOP > KBOT) return;
            // *     ... nor for an empty deflation window. ====
            if (NW < 1) return;
            // *
            // *     ==== Machine constants ====
            // *
            SAFMIN = _dlamch.Run("SAFE MINIMUM");
            SAFMAX = ONE / SAFMIN;
            _dlabad.Run(ref SAFMIN, ref SAFMAX);
            ULP = _dlamch.Run("PRECISION");
            SMLNUM = SAFMIN * (Convert.ToDouble(N) / ULP);
            // *
            // *     ==== Setup deflation window ====
            // *
            JW = Math.Min(NW, KBOT - KTOP + 1);
            KWTOP = KBOT - JW + 1;
            if (KWTOP == KTOP)
            {
                S = ZERO;
            }
            else
            {
                S = H[KWTOP+(KWTOP - 1) * LDH + o_h];
            }
            // *
            if (KBOT == KWTOP)
            {
                // *
                // *        ==== 1-by-1 deflation window: not much to do ====
                // *
                SR[KWTOP + o_sr] = H[KWTOP+KWTOP * LDH + o_h];
                SI[KWTOP + o_si] = ZERO;
                NS = 1;
                ND = 0;
                if (Math.Abs(S) <= Math.Max(SMLNUM, ULP * Math.Abs(H[KWTOP+KWTOP * LDH + o_h])))
                {
                    NS = 0;
                    ND = 1;
                    if (KWTOP > KTOP) H[KWTOP+(KWTOP - 1) * LDH + o_h] = ZERO;
                }
                return;
            }
            // *
            // *     ==== Convert to spike-triangular form.  (In case of a
            // *     .    rare QR failure, this routine continues to do
            // *     .    aggressive early deflation using that part of
            // *     .    the deflation window that converged using INFQR
            // *     .    here and there to keep track.) ====
            // *
            _dlacpy.Run("U", JW, JW, H, KWTOP+KWTOP * LDH + o_h, LDH, ref T, offset_t
                             , LDT);
            _dcopy.Run(JW - 1, H, KWTOP + 1+KWTOP * LDH + o_h, LDH + 1, ref T, 2+1 * LDT + o_t, LDT + 1);
            // *
            _dlaset.Run("A", JW, JW, ZERO, ONE, ref V, offset_v
                             , LDV);
            _dlahqr.Run(true, true, JW, 1, JW, ref T, offset_t
                             , LDT, ref SR, KWTOP + o_sr, ref SI, KWTOP + o_si, 1, JW, ref V, offset_v
                             , LDV, ref INFQR);
            // *
            // *     ==== DTREXC needs a clean margin near the diagonal ====
            // *
            for (J = 1; J <= JW - 3; J++)
            {
                T[J + 2+J * LDT + o_t] = ZERO;
                T[J + 3+J * LDT + o_t] = ZERO;
            }
            if (JW > 2) T[JW+(JW - 2) * LDT + o_t] = ZERO;
            // *
            // *     ==== Deflation detection loop ====
            // *
            NS = JW;
            ILST = INFQR + 1;
        LABEL20:;
            if (ILST <= NS)
            {
                if (NS == 1)
                {
                    BULGE = false;
                }
                else
                {
                    BULGE = T[NS+(NS - 1) * LDT + o_t] != ZERO;
                }
                // *
                // *        ==== Small spike tip test for deflation ====
                // *
                if (!BULGE)
                {
                    // *
                    // *           ==== Real eigenvalue ====
                    // *
                    FOO = Math.Abs(T[NS+NS * LDT + o_t]);
                    if (FOO == ZERO) FOO = Math.Abs(S);
                    if (Math.Abs(S * V[1+NS * LDV + o_v]) <= Math.Max(SMLNUM, ULP * FOO))
                    {
                        // *
                        // *              ==== Deflatable ====
                        // *
                        NS -= 1;
                    }
                    else
                    {
                        // *
                        // *              ==== Undeflatable.   Move it up out of the way.
                        // *              .    (DTREXC can not fail in this case.) ====
                        // *
                        IFST = NS;
                        _dtrexc.Run("V", JW, ref T, offset_t, LDT, ref V, offset_v, LDV
                                         , ref IFST, ref ILST, ref WORK, offset_work, ref INFO);
                        ILST += 1;
                    }
                }
                else
                {
                    // *
                    // *           ==== Complex conjugate pair ====
                    // *
                    FOO = Math.Abs(T[NS+NS * LDT + o_t]) + Math.Sqrt(Math.Abs(T[NS+(NS - 1) * LDT + o_t])) * Math.Sqrt(Math.Abs(T[NS - 1+NS * LDT + o_t]));
                    if (FOO == ZERO) FOO = Math.Abs(S);
                    if (Math.Max(Math.Abs(S * V[1+NS * LDV + o_v]), Math.Abs(S * V[1+(NS - 1) * LDV + o_v])) <= Math.Max(SMLNUM, ULP * FOO))
                    {
                        // *
                        // *              ==== Deflatable ====
                        // *
                        NS -= 2;
                    }
                    else
                    {
                        // *
                        // *              ==== Undflatable. Move them up out of the way.
                        // *              .    Fortunately, DTREXC does the right thing with
                        // *              .    ILST in case of a rare exchange failure. ====
                        // *
                        IFST = NS;
                        _dtrexc.Run("V", JW, ref T, offset_t, LDT, ref V, offset_v, LDV
                                         , ref IFST, ref ILST, ref WORK, offset_work, ref INFO);
                        ILST += 2;
                    }
                }
                // *
                // *        ==== End deflation detection loop ====
                // *
                goto LABEL20;
            }
            // *
            // *        ==== Return to Hessenberg form ====
            // *
            if (NS == 0) S = ZERO;
            // *
            if (NS < JW)
            {
                // *
                // *        ==== sorting diagonal blocks of T improves accuracy for
                // *        .    graded matrices.  Bubble sort deals well with
                // *        .    exchange failures. ====
                // *
                SORTED = false;
                I = NS + 1;
            LABEL30:;
                if (SORTED) goto LABEL50;
                SORTED = true;
                // *
                KEND = I - 1;
                I = INFQR + 1;
                if (I == NS)
                {
                    K = I + 1;
                }
                else
                {
                    if (T[I + 1+I * LDT + o_t] == ZERO)
                    {
                        K = I + 1;
                    }
                    else
                    {
                        K = I + 2;
                    }
                }
            LABEL40:;
                if (K <= KEND)
                {
                    if (K == I + 1)
                    {
                        EVI = Math.Abs(T[I+I * LDT + o_t]);
                    }
                    else
                    {
                        EVI = Math.Abs(T[I+I * LDT + o_t]) + Math.Sqrt(Math.Abs(T[I + 1+I * LDT + o_t])) * Math.Sqrt(Math.Abs(T[I+(I + 1) * LDT + o_t]));
                    }
                    // *
                    if (K == KEND)
                    {
                        EVK = Math.Abs(T[K+K * LDT + o_t]);
                    }
                    else
                    {
                        if (T[K + 1+K * LDT + o_t] == ZERO)
                        {
                            EVK = Math.Abs(T[K+K * LDT + o_t]);
                        }
                        else
                        {
                            EVK = Math.Abs(T[K+K * LDT + o_t]) + Math.Sqrt(Math.Abs(T[K + 1+K * LDT + o_t])) * Math.Sqrt(Math.Abs(T[K+(K + 1) * LDT + o_t]));
                        }
                    }
                    // *
                    if (EVI >= EVK)
                    {
                        I = K;
                    }
                    else
                    {
                        SORTED = false;
                        IFST = I;
                        ILST = K;
                        _dtrexc.Run("V", JW, ref T, offset_t, LDT, ref V, offset_v, LDV
                                         , ref IFST, ref ILST, ref WORK, offset_work, ref INFO);
                        if (INFO == 0)
                        {
                            I = ILST;
                        }
                        else
                        {
                            I = K;
                        }
                    }
                    if (I == KEND)
                    {
                        K = I + 1;
                    }
                    else
                    {
                        if (T[I + 1+I * LDT + o_t] == ZERO)
                        {
                            K = I + 1;
                        }
                        else
                        {
                            K = I + 2;
                        }
                    }
                    goto LABEL40;
                }
                goto LABEL30;
            LABEL50:;
            }
            // *
            // *     ==== Restore shift/eigenvalue array from T ====
            // *
            I = JW;
        LABEL60:;
            if (I >= INFQR + 1)
            {
                if (I == INFQR + 1)
                {
                    SR[KWTOP + I - 1 + o_sr] = T[I+I * LDT + o_t];
                    SI[KWTOP + I - 1 + o_si] = ZERO;
                    I -= 1;
                }
                else
                {
                    if (T[I+(I - 1) * LDT + o_t] == ZERO)
                    {
                        SR[KWTOP + I - 1 + o_sr] = T[I+I * LDT + o_t];
                        SI[KWTOP + I - 1 + o_si] = ZERO;
                        I -= 1;
                    }
                    else
                    {
                        AA = T[I - 1+(I - 1) * LDT + o_t];
                        CC = T[I+(I - 1) * LDT + o_t];
                        BB = T[I - 1+I * LDT + o_t];
                        DD = T[I+I * LDT + o_t];
                        _dlanv2.Run(ref AA, ref BB, ref CC, ref DD, ref SR[KWTOP + I - 2 + o_sr], ref SI[KWTOP + I - 2 + o_si]
                                         , ref SR[KWTOP + I - 1 + o_sr], ref SI[KWTOP + I - 1 + o_si], ref CS, ref SN);
                        I -= 2;
                    }
                }
                goto LABEL60;
            }
            // *
            if (NS < JW || S == ZERO)
            {
                if (NS > 1 && S != ZERO)
                {
                    // *
                    // *           ==== Reflect spike back into lower triangle ====
                    // *
                    _dcopy.Run(NS, V, offset_v, LDV, ref WORK, offset_work, 1);
                    BETA = WORK[1 + o_work];
                    _dlarfg.Run(NS, ref BETA, ref WORK, 2 + o_work, 1, ref TAU);
                    WORK[1 + o_work] = ONE;
                    // *
                    _dlaset.Run("L", JW - 2, JW - 2, ZERO, ZERO, ref T, 3+1 * LDT + o_t
                                     , LDT);
                    // *
                    _dlarf.Run("L", NS, JW, WORK, offset_work, 1, TAU
                                    , ref T, offset_t, LDT, ref WORK, JW + 1 + o_work);
                    _dlarf.Run("R", NS, NS, WORK, offset_work, 1, TAU
                                    , ref T, offset_t, LDT, ref WORK, JW + 1 + o_work);
                    _dlarf.Run("R", JW, NS, WORK, offset_work, 1, TAU
                                    , ref V, offset_v, LDV, ref WORK, JW + 1 + o_work);
                    // *
                    _dgehrd.Run(JW, 1, NS, ref T, offset_t, LDT, ref WORK, offset_work
                                     , ref WORK, JW + 1 + o_work, LWORK - JW, ref INFO);
                }
                // *
                // *        ==== Copy updated reduced window into place ====
                // *
                if (KWTOP > 1) H[KWTOP+(KWTOP - 1) * LDH + o_h] = S * V[1+1 * LDV + o_v];
                _dlacpy.Run("U", JW, JW, T, offset_t, LDT, ref H, KWTOP+KWTOP * LDH + o_h
                                 , LDH);
                _dcopy.Run(JW - 1, T, 2+1 * LDT + o_t, LDT + 1, ref H, KWTOP + 1+KWTOP * LDH + o_h, LDH + 1);
                // *
                // *        ==== Accumulate orthogonal matrix in order update
                // *        .    H and Z, if requested.  (A modified version
                // *        .    of  DORGHR that accumulates block Householder
                // *        .    transformations into V directly might be
                // *        .    marginally more efficient than the following.) ====
                // *
                if (NS > 1 && S != ZERO)
                {
                    _dorghr.Run(JW, 1, NS, ref T, offset_t, LDT, WORK, offset_work
                                     , ref WORK, JW + 1 + o_work, LWORK - JW, ref INFO);
                    _dgemm.Run("N", "N", JW, NS, NS, ONE
                                    , V, offset_v, LDV, T, offset_t, LDT, ZERO, ref WV, offset_wv
                                    , LDWV);
                    _dlacpy.Run("A", JW, NS, WV, offset_wv, LDWV, ref V, offset_v
                                     , LDV);
                }
                // *
                // *        ==== Update vertical slab in H ====
                // *
                if (WANTT)
                {
                    LTOP = 1;
                }
                else
                {
                    LTOP = KTOP;
                }
                for (KROW = LTOP; NV >= 0 ? KROW <= KWTOP - 1 : KROW >= KWTOP - 1; KROW += NV)
                {
                    KLN = Math.Min(NV, KWTOP - KROW);
                    _dgemm.Run("N", "N", KLN, JW, JW, ONE
                                    , H, KROW+KWTOP * LDH + o_h, LDH, V, offset_v, LDV, ZERO, ref WV, offset_wv
                                    , LDWV);
                    _dlacpy.Run("A", KLN, JW, WV, offset_wv, LDWV, ref H, KROW+KWTOP * LDH + o_h
                                     , LDH);
                }
                // *
                // *        ==== Update horizontal slab in H ====
                // *
                if (WANTT)
                {
                    for (KCOL = KBOT + 1; NH >= 0 ? KCOL <= N : KCOL >= N; KCOL += NH)
                    {
                        KLN = Math.Min(NH, N - KCOL + 1);
                        _dgemm.Run("C", "N", JW, KLN, JW, ONE
                                        , V, offset_v, LDV, H, KWTOP+KCOL * LDH + o_h, LDH, ZERO, ref T, offset_t
                                        , LDT);
                        _dlacpy.Run("A", JW, KLN, T, offset_t, LDT, ref H, KWTOP+KCOL * LDH + o_h
                                         , LDH);
                    }
                }
                // *
                // *        ==== Update vertical slab in Z ====
                // *
                if (WANTZ)
                {
                    for (KROW = ILOZ; NV >= 0 ? KROW <= IHIZ : KROW >= IHIZ; KROW += NV)
                    {
                        KLN = Math.Min(NV, IHIZ - KROW + 1);
                        _dgemm.Run("N", "N", KLN, JW, JW, ONE
                                        , Z, KROW+KWTOP * LDZ + o_z, LDZ, V, offset_v, LDV, ZERO, ref WV, offset_wv
                                        , LDWV);
                        _dlacpy.Run("A", KLN, JW, WV, offset_wv, LDWV, ref Z, KROW+KWTOP * LDZ + o_z
                                         , LDZ);
                    }
                }
            }
            // *
            // *     ==== Return the number of deflations ... ====
            // *
            ND = JW - NS;
            // *
            // *     ==== ... and the number of shifts. (Subtracting
            // *     .    INFQR from the spike length takes care
            // *     .    of the case of a rare QR failure while
            // *     .    calculating eigenvalues of the deflation
            // *     .    window.)  ====
            // *
            NS -= INFQR;
            // *
            // *      ==== Return optimal workspace. ====
            // *
            WORK[1 + o_work] = Convert.ToDouble(LWKOPT);
            // *
            // *     ==== End of DLAQR2 ====
            // *

            #endregion

        }
    }
}
