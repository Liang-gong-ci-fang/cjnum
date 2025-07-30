package cjnum.lapack

import cjnum.blas.*

public interface LFloat64 {
    func dgecon(norm: MatrixNorm, n: Int64, a: Array<Float64>, lda: Int64, anorm: Float64, work: Array<Float64>,
        iwork: Array<Int64>): Float64
    func dgeev(jobvl: LeftEVJob, jobvr: RightEVJob, n: Int64, a: Array<Float64>, lda: Int64, wr: Array<Float64>,
        wi: Array<Float64>, vl: Array<Float64>, ldvl: Int64, vr: Array<Float64>, ldvr: Int64, work: Array<Float64>,
        lwork: Int64): Int64
    func dgels(trans: Transpose, m: Int64, n: Int64, nrhs: Int64, a: Array<Float64>, lda: Int64, b: Array<Float64>,
        ldb: Int64, work: Array<Float64>, lwork: Int64): Bool
    func dgelqf(m: Int64, n: Int64, a: Array<Float64>, lda: Int64, tau: Array<Float64>, work: Array<Float64>,
        lwork: Int64): Unit
    func dgeqp3(m: Int64, n: Int64, a: Array<Float64>, lda: Int64, jpvt: Array<Int64>, tau: Array<Float64>,
        work: Array<Float64>, lwork: Int64): Unit
    func dgeqrf(m: Int64, n: Int64, a: Array<Float64>, lda: Int64, tau: Array<Float64>, work: Array<Float64>,
        lwork: Int64): Unit
    func dgesvd(jobU: SVDJob, jobVT: SVDJob, m: Int64, n: Int64, a: Array<Float64>, lda: Int64, s: Array<Float64>,
        u: Array<Float64>, ldu: Int64, vt: Array<Float64>, ldvt: Int64, work: Array<Float64>, lwork: Int64): Bool
    func dgetrf(m: Int64, n: Int64, a: Array<Float64>, lda: Int64, ipiv: Array<Int64>): Bool
    func dgetri(n: Int64, a: Array<Float64>, lda: Int64, ipiv: Array<Int64>, work: Array<Float64>, lwork: Int64): Bool
    func dgetrs(trans: Transpose, n: Int64, nrhs: Int64, a: Array<Float64>, lda: Int64, ipiv: Array<Int64>,
        b: Array<Float64>, ldb: Int64): Unit
    func dggsvd3(jobU: GSVDJob, jobV: GSVDJob, jobQ: GSVDJob, m: Int64, n: Int64, p: Int64, a: Array<Float64>,
        lda: Int64, b: Array<Float64>, ldb: Int64, alpha: Array<Float64>, beta: Array<Float64>, u: Array<Float64>,
        ldu: Int64, v: Array<Float64>, ldv: Int64, q: Array<Float64>, ldq: Int64, work: Array<Float64>, lwork: Int64,
        iwork: Array<Int64>): (Int64, Int64, Bool)
    func dlantr(norm: MatrixNorm, uplo: Uplo, diag: Diag, m: Int64, n: Int64, a: Array<Float64>, lda: Int64,
        work: Array<Float64>): Float64
    func dlange(norm: MatrixNorm, m: Int64, n: Int64, a: Array<Float64>, lda: Int64, work: Array<Float64>): Float64
    func dlansy(norm: MatrixNorm, uplo: Uplo, n: Int64, a: Array<Float64>, lda: Int64, work: Array<Float64>): Float64
    func dlapmr(forward: Bool, m: Int64, n: Int64, x: Array<Float64>, ldx: Int64, k: Array<Int64>): Unit
    func dlapmt(forward: Bool, m: Int64, n: Int64, x: Array<Float64>, ldx: Int64, k: Array<Int64>): Unit
    func dorgqr(m: Int64, n: Int64, k: Int64, a: Array<Float64>, lda: Int64, tau: Array<Float64>, work: Array<Float64>,
        lwork: Int64): Unit
    func dormqr(side: Side, trans: Transpose, m: Int64, n: Int64, k: Int64, a: Array<Float64>, lda: Int64,
        tau: Array<Float64>, c: Array<Float64>, ldc: Int64, work: Array<Float64>, lwork: Int64): Unit
    func dorglq(m: Int64, n: Int64, k: Int64, a: Array<Float64>, lda: Int64, tau: Int64, work: Array<Float64>,
        lwork: Int64): Unit
    func dormlq(side: Side, trans: Transpose, m: Int64, n: Int64, k: Int64, a: Array<Float64>, lda: Int64,
        tau: Array<Float64>, c: Array<Float64>, ldc: Int64, work: Array<Float64>, lwork: Int64): Unit
    func dpbcon(uplo: Uplo, n: Int64, kd: Int64, ab: Array<Float64>, ldab: Int64, anorm: Float64, work: Array<Float64>,
        iwork: Array<Int64>): Float64
    func dpbtrf(uplo: Uplo, n: Int64, kd: Int64, ab: Array<Float64>, ldab: Int64): Bool
    func dpbtrs(uplo: Uplo, n: Int64, kd: Int64, nrhs: Int64, ab: Array<Float64>, ldab: Int64, b: Array<Float64>,
        ldb: Int64): Unit
    func dpocon(uplo: Uplo, n: Int64, a: Array<Float64>, lda: Int64, anorm: Float64, work: Array<Float64>,
        iwork: Array<Int64>): Float64
    func dpotrf(ul: Uplo, n: Int64, a: Array<Float64>, lda: Int64): Bool
    func dpotri(ul: Uplo, n: Int64, a: Array<Float64>, lda: Int64): Bool
    func dpotrs(ul: Uplo, n: Int64, nrhs: Int64, a: Array<Float64>, lda: Int64, b: Array<Float64>, ldb: Int64): Unit
    func dpstrf(uplo: Uplo, n: Int64, a: Array<Float64>, lda: Int64, piv: Array<Int64>, tol: Float64,
        work: Array<Float64>): (Int64, Bool)
    func dsyev(jobz: EVJob, uplo: Uplo, n: Int64, a: Array<Float64>, lda: Int64, w: Array<Float64>,
        work: Array<Float64>, lwork: Int64): Bool
    func dtbtrs(uplo: Uplo, trans: Transpose, diag: Diag, n: Int64, kd: Int64, nrhs: Int64, a: Array<Float64>,
        lda: Int64, b: Array<Float64>, ldb: Int64): Bool
    func dtrcon(norm: MatrixNorm, uplo: Uplo, diag: Diag, n: Int64, a: Array<Float64>, lda: Int64, work: Array<Float64>,
        iwork: Array<Int64>): Float64
    func dtrtri(uplo: Uplo, diag: Diag, n: Int64, a: Array<Float64>, lda: Int64): Bool
    func dtrtrs(uplo: Uplo, trans: Transpose, diag: Diag, n: Int64, nrhs: Int64, a: Array<Float64>, lda: Int64,
        b: Array<Float64>, ldb: Int64): Bool
    func drscl(n: Int64, a: Float64, x: Array<Float64>, incX: Int64): Unit
    func dlassq(n: Int64, x: Array<Float64>, incX: Int64, scale: Float64, sumsq: Float64): (Float64, Float64)
    func dlatrs(uplo: Uplo, trans: Transpose, diag: Diag, normin: Bool, n: Int64, a: Array<Float64>, lda: Int64,
        x: Array<Float64>, cnorm: Array<Float64>): Float64
    func dlacn2(n: Int64, v: Array<Float64>, x: Array<Float64>, isgn: Array<Int64>, est: Float64, kase: Int64,
        isave: Array<Int64>): (Float64, Int64)
    func dtrti2(uplo: Uplo, diag: Diag, n: Int64, a: Array<Float64>, lda: Int64): Unit
    func dgetf2(m: Int64, n: Int64, a: Array<Float64>, lda: Int64, ipiv: Array<Int64>): Bool
    func dlaswp(n: Int64, a: Array<Float64>, lda: Int64, k1: Int64, k2: Int64, ipiv: Array<Int64>, incX: Int64): Unit
    func dlartg(f: Float64, g: Float64): (Float64, Float64, Float64)
    func dlas2(f: Float64, g: Float64, h: Float64): (Float64, Float64)
    func dlae2(a: Float64, b: Float64, c: Float64): (Float64, Float64)
    func dlacpy(uplo: Uplo, m: Int64, n: Int64, a: Array<Float64>, lda: Int64, b: Array<Float64>, ldb: Int64): Unit
    func iladlc(m: Int64, n: Int64, a: Array<Float64>, lda: Int64): Int64
    func iladlr(m: Int64, n: Int64, a: Array<Float64>, lda: Int64): Int64
    func dlaset(uplo: Uplo, m: Int64, n: Int64, alpha: Float64, beta: Float64, a: Array<Float64>, lda: Int64): Unit
    func dlaqr1(n: Int64, h: Array<Float64>, ldh: Int64, sr1: Float64, si1: Float64, sr2: Float64, si2: Float64,
        v: Array<Float64>): Unit
    func dptts2(n: Int64, nrhs: Int64, d: Array<Float64>, e: Array<Float64>, b: Array<Float64>, ldb: Int64): Unit
    func dpttrs(n: Int64, nrhs: Int64, d: Array<Float64>, e: Array<Float64>, b: Array<Float64>, ldb: Int64): Unit
    func dpttrf(n: Int64, d: Array<Float64>, e: Array<Float64>): Bool
    func dlarfg(n: Int64, alpha: Float64, x: Array<Float64>, incX: Int64): (Float64, Float64)
    func dlarf(side: Side, m: Int64, n: Int64, v: Array<Float64>, incv: Int64, tau: Float64, c: Array<Float64>,
        ldc: Int64, work: Array<Float64>): Unit
    func dlanhs(norm: MatrixNorm, n: Int64, a: Array<Float64>, lda: Int64, work: Array<Float64>): Float64
    func dgeql2(m: Int64, n: Int64, a: Array<Float64>, lda: Int64, tau: Array<Float64>, work: Array<Float64>): Unit
    func dlanv2(a: Float64, b: Float64, c: Float64, d: Float64): (Float64, Float64, Float64, Float64, Float64, Float64, 
        Float64, Float64, Float64, Float64)
    func dgebd2(m: Int64, n: Int64, a: Array<Float64>, lda: Int64, d: Array<Float64>, e: Array<Float64>,
        tauQ: Array<Float64>, tauP: Array<Float64>, work: Array<Float64>): Unit
    func dlarfx(side: Side, m: Int64, n: Int64, v: Array<Float64>, tau: Float64, c: Array<Float64>, ldc: Int64,
        work: Array<Float64>): Unit
    func dlansb(norm: MatrixNorm, uplo: Uplo, n: Int64, kd: Int64, ab: Array<Float64>, ldab: Int64, work: Array<Float64>): Float64
    func dlanst(norm: MatrixNorm, n: Int64, d: Array<Float64>, e: Array<Float64>): Float64
    func dlantb(norm: MatrixNorm, uplo: Uplo, diag: Diag, n: Int64, k: Int64, a: Array<Float64>, lda: Int64,
        work: Array<Float64>): Float64
    func dlasr(side: Side, pivot: Pivot, direct: Direct, m: Int64, n: Int64, c: Array<Float64>, s: Array<Float64>,
        a: Array<Float64>, lda: Int64): Unit
    func dlasy2(tranl: Bool, tranr: Bool, isgn: Int64, n1: Int64, n2: Int64, tl: Array<Float64>, ldtl: Int64,
        tr: Array<Float64>, ldtr: Int64, b: Array<Float64>, ldb: Int64, x: Array<Float64>, ldx: Int64): (Float64, 
        Float64, Bool)
    func dlascl(kind: MatrixType, kl: Int64, ku: Int64, cfrom: Float64, cto: Float64, m: Int64, n: Int64,
        a: Array<Float64>, lda: Int64): Unit
    func dlaexc(wantq: Bool, n: Int64, t: Array<Float64>, ldt: Int64, q: Array<Float64>, ldq: Int64, j1: Int64,
        n1: Int64, n2: Int64, work: Array<Float64>): Bool
    func dlasv2(f: Float64, g: Float64, h: Float64): (ssmin: Float64, ssmax: Float64, snr: Float64, csr: Float64, 
        snl: Float64, csl: Float64)
    func dlauu2(uplo: Uplo, n: Int64, a: Array<Float64>, lda: Int64): Unit
    func dpotf2(uplo: Uplo, n: Int64, a: Array<Float64>, lda: Int64): Bool
    func dpbtf2(uplo: Uplo, n: Int64, kd: Int64, ab: Array<Float64>, ldab: Int64): Bool
    func dgetc2(n: Int64, a: Array<Float64>, lda: Int64, ipiv: Array<Int64>, jpiv: Array<Int64>): Int64
    func dgesc2(n: Int64, a: Array<Float64>, lda: Int64, rhs: Array<Float64>, ipiv: Array<Int64>, jpiv: Array<Int64>): Float64
    func dgghrd(compq: OrthoComp, compz: OrthoComp, n: Int64, ilo: Int64, ihi: Int64, a: Array<Float64>, lda: Int64,
        b: Array<Float64>, ldb: Int64, q: Array<Float64>, ldq: Int64, z: Array<Float64>, ldz: Int64): Unit
    func dgtsv(n: Int64, nrhs: Int64, dl: Array<Float64>, d: Array<Float64>, du: Array<Float64>, b: Array<Float64>,
        ldb: Int64): Bool
    func dgehd2(n: Int64, ilo: Int64, ihi: Int64, a: Array<Float64>, lda: Int64, tau: Array<Float64>,
        work: Array<Float64>): Unit
    func dlatrd(uplo: Uplo, n: Int64, nb: Int64, a: Array<Float64>, lda: Int64, e: Array<Float64>, tau: Array<Float64>,
        w: Array<Float64>, ldw: Int64): Unit
    func dlahqr(wantt: Bool, wantz: Bool, n: Int64, ilo: Int64, ihi: Int64, h: Array<Float64>, ldh: Int64,
        wr: Array<Float64>, wi: Array<Float64>, iloz: Int64, ihiz: Int64, z: Array<Float64>, ldz: Int64): Int64
    func dgeqr2(m: Int64, n: Int64, a: Array<Float64>, lda: Int64, tau: Array<Float64>, work: Array<Float64>): Unit
    func dlaqr5(wantt: Bool, wantz: Bool, kacc22: Int64, n: Int64, ktop: Int64, kbot: Int64, nshfts: Int64,
        sr: Array<Float64>, si: Array<Float64>, h: Array<Float64>, ldh: Int64, iloz: Int64, ihiz: Int64,
        z: Array<Float64>, ldz: Int64, v: Array<Float64>, ldv: Int64, u: Array<Float64>, ldu: Int64, nv: Int64,
        wv: Array<Float64>, ldwv: Int64, nh: Int64, wh: Array<Float64>, ldwh: Int64): Unit
    func dlarft(direct: Direct, store: StoreV, n: Int64, k: Int64, v: Array<Float64>, ldv: Int64, tau: Array<Float64>,
        t: Array<Float64>, ldt: Int64): Unit
    func dorg2l(m: Int64, n: Int64, k: Int64, a: Array<Float64>, lda: Int64, tau: Array<Float64>, work: Array<Float64>): Unit
    func dlarfb(side: Side, trans: Transpose, direct: Direct, store: StoreV, m: Int64, n: Int64, k: Int64,
        v: Array<Float64>, ldv: Int64, t: Array<Float64>, ldt: Int64, c: Array<Float64>, ldc: Int64,
        work: Array<Float64>, ldwork: Int64): Unit
    func dorgql(m: Int64, n: Int64, k: Int64, a: Array<Float64>, lda: Int64, tau: Array<Float64>, work: Array<Float64>,
        lwork: Int64): Unit
    func dlaqr04(wantt: Bool, wantz: Bool, n: Int64, ilo: Int64, ihi: Int64, h: Array<Float64>, ldh: Int64,
        wr: Array<Float64>, wi: Array<Float64>, iloz: Int64, ihiz: Int64, z: Array<Float64>, ldz: Int64,
        work: Array<Float64>, lwork: Int64, recur: Int64): Int64
    func dlaqr23(wantt: Bool, wantz: Bool, n: Int64, ktop: Int64, kbot: Int64, nw: Int64, h: Array<Float64>, ldh: Int64,
        iloz: Int64, ihiz: Int64, z: Array<Float64>, ldz: Int64, sr: Array<Float64>, si: Array<Float64>,
        v: Array<Float64>, ldv: Int64, nh: Int64, t: Array<Float64>, ldt: Int64, nv: Int64, wv: Array<Float64>,
        ldwv: Int64, work: Array<Float64>, lwork: Int64, recur: Int64): (Int64, Int64)
    func dtrexc(compq: UpdateSchurComp, n: Int64, t: Array<Float64>, ldt: Int64, q: Array<Float64>, ldq: Int64,
        ifst: Int64, ilst: Int64, work: Array<Float64>): (Int64, Int64, Bool)
    func dlasrt(s: Sort, n: Int64, d: Array<Float64>): Unit
    func dlaev2(a: Float64, b: Float64, c: Float64): (Float64, Float64, Float64, Float64)
    func dpstf2(uplo: Uplo, n: Int64, a: Array<Float64>, lda: Int64, piv: Array<Int64>, tol: Float64,
        work: Array<Float64>): (Int64, Bool)
    func dsytd2(uplo: Uplo, n: Int64, a: Array<Float64>, lda: Int64, d: Array<Float64>, e: Array<Float64>,
        tau: Array<Float64>): Unit
    func dlahr2(n: Int64, k: Int64, nb: Int64, a: Array<Float64>, lda: Int64, tau: Array<Float64>, t: Array<Float64>,
        ldt: Int64, y: Array<Float64>, ldy: Int64): Unit
    func dlag2(a: Array<Float64>, lda: Int64, b: Array<Float64>, ldb: Int64): (Float64, Float64, Float64, Float64, 
        Float64)
    func dlags2(upper: Bool, a1: Float64, a2: Float64, a3: Float64, b1: Float64, b2: Float64, b3: Float64): (Float64, 
        Float64, Float64, Float64, Float64, Float64)
    func dlagtm(trans: Transpose, m: Int64, n: Int64, alpha: Float64, dl: Array<Float64>, d: Array<Float64>,
        du: Array<Float64>, b: Array<Float64>, ldb: Int64, beta: Float64, c: Array<Float64>, ldc: Int64): Unit
    func dlangt(norm: MatrixNorm, n: Int64, dl: Array<Float64>, d: Array<Float64>, du: Array<Float64>): Float64
    func dlaln2(trans: Bool, na: Int64, nw: Int64, smin: Float64, ca: Float64, a: Array<Float64>, lda: Int64,
        d1: Float64, d2: Float64, b: Array<Float64>, ldb: Int64, wr: Float64, wi: Float64, x: Array<Float64>, ldx: Int64): (Float64, 
        Float64, Bool)
    func dgebak(job: BalanceJob, side: EVSide, n: Int64, ilo: Int64, ihi: Int64, scale: Array<Float64>, m: Int64,
        v: Array<Float64>, ldv: Int64): Unit
    func dgelq2(m: Int64, n: Int64, a: Array<Float64>, lda: Int64, tau: Array<Float64>, work: Array<Float64>): Unit
    func dgerq2(m: Int64, n: Int64, a: Array<Float64>, lda: Int64, tau: Array<Float64>, work: Array<Float64>): Unit
    func dlangb(norm: MatrixNorm, m: Int64, n: Int64, kl: Int64, ku: Int64, ab: Array<Float64>, ldab: Int64): Float64
    func dlatbs(uplo: Uplo, trans: Transpose, diag: Diag, normin: Bool, n: Int64, kd: Int64, ab: Array<Float64>,
        ldab: Int64, x: Array<Float64>, cnorm: Array<Float64>): Float64
    func dgebrd(m: Int64, n: Int64, a: Array<Float64>, lda: Int64, d: Array<Float64>, e: Array<Float64>,
        tauQ: Array<Float64>, tauP: Array<Float64>, work: Array<Float64>, lwork: Int64): Unit
    func dlabrd(m: Int64, n: Int64, nb: Int64, a: Array<Float64>, lda: Int64, d: Array<Float64>, e: Array<Float64>,
        tauQ: Array<Float64>, tauP: Array<Float64>, x: Array<Float64>, ldx: Int64, y: Array<Float64>, ldy: Int64): Unit
    func dlasq1(n: Int64, d: Array<Float64>, e: Array<Float64>, work: Array<Float64>): Int64
    func dlasq2(n: Int64, z: Array<Float64>): Int64
    func dorgtr(uplo: Uplo, n: Int64, a: Array<Float64>, lda: Int64, tau: Array<Float64>, work: Array<Float64>,
        lwork: Int64): Unit
    func dsteqr(compz: EVComp, n: Int64, d: Array<Float64>, e: Array<Float64>, z: Array<Float64>, ldz: Int64,
        work: Array<Float64>): Bool
    func dsterf(n: Int64, d: Array<Float64>, e: Array<Float64>): Bool
    func dsytrd(uplo: Uplo, n: Int64, a: Array<Float64>, lda: Int64, d: Array<Float64>, e: Array<Float64>,
        tau: Array<Float64>, work: Array<Float64>, lwork: Int64): Unit
    func dormhr(side: Side, trans: Transpose, m: Int64, n: Int64, ilo: Int64, ihi: Int64, a: Array<Float64>, lda: Int64,
        tau: Array<Float64>, c: Array<Float64>, ldc: Int64, work: Array<Float64>, lwork: Int64): Unit
    func dgehrd(n: Int64, ilo: Int64, ihi: Int64, a: Array<Float64>, lda: Int64, tau: Array<Float64>,
        work: Array<Float64>, lwork: Int64): Unit
    func dorg2r(m: Int64, n: Int64, k: Int64, a: Array<Float64>, lda: Int64, tau: Array<Float64>, work: Array<Float64>): Unit
    func dorm2r(side: Side, trans: Transpose, m: Int64, n: Int64, k: Int64, a: Array<Float64>, lda: Int64,
        tau: Array<Float64>, c: Array<Float64>, ldc: Int64, work: Array<Float64>): Unit
    func dtrevc3(side: EVSide, howmny: EVHowMany, selected: Array<Bool>, n: Int64, t: Array<Float64>, ldt: Int64,
        vl: Array<Float64>, ldvl: Int64, vr: Array<Float64>, ldvr: Int64, mm: Int64, work: Array<Float64>, lwork: Int64): Int64
    func dormr2(side: Side, trans: Transpose, m: Int64, n: Int64, k: Int64, a: Array<Float64>, lda: Int64,
        tau: Array<Float64>, c: Array<Float64>, ldc: Int64, work: Array<Float64>): Unit
    func dgerqf(m: Int64, n: Int64, a: Array<Float64>, lda: Int64, tau: Array<Float64>, work: Array<Float64>,
        lwork: Int64): Unit
    func dorml2(side: Side, trans: Transpose, m: Int64, n: Int64, k: Int64, a: Array<Float64>, lda: Int64,
        tau: Array<Float64>, c: Array<Float64>, ldc: Int64, work: Array<Float64>): Unit
    func dormbr(vect: ApplyOrtho, side: Side, trans: Transpose, m: Int64, n: Int64, k: Int64, a: Array<Float64>,
        lda: Int64, tau: Array<Float64>, c: Array<Float64>, ldc: Int64, work: Array<Float64>, lwork: Int64): Unit
}

public type Direct = UInt8

public const Forward: Direct = 'F'
public const Backward: Direct = 'B'

public type Sort = UInt8

public const SortIncreasing: Sort = 'I'
public const SortDecreasing: Sort = 'D'

public type StoreV = UInt8

public const ColumnWise: StoreV = 'C'
public const RowWise: StoreV = 'R'

public type MatrixNorm = UInt8

public const MaxAbs: MatrixNorm = 'M'
public const MaxColumnSum: MatrixNorm = 'O'
public const MaxRowSum: MatrixNorm = 'I'
public const Frobenius: MatrixNorm = 'F'

public type MatrixType = UInt8

public const General: MatrixType = 'G'
public const UpperTri: MatrixType = 'U'
public const LowerTri: MatrixType = 'L'

public type Pivot = UInt8

public const Variable: Pivot = 'V'
public const Top: Pivot = 'T'
public const Bottom: Pivot = 'B'

public type ApplyOrtho = UInt8

public const ApplyP: ApplyOrtho = 'P'
public const ApplyQ: ApplyOrtho = 'Q'

public type GenOrtho = UInt8

public const GeneratePT: GenOrtho = 'P'
public const GenerateQ: GenOrtho = 'Q'

public type SVDJob = UInt8

public const SVDAll: SVDJob = 'A'
public const SVDStore: SVDJob = 'S'
public const SVDOverwrite: SVDJob = 'O'
public const SVDNone: SVDJob = 'N'

public type GSVDJob = UInt8

public const GSVDU: GSVDJob = 'U'
public const GSVDV: GSVDJob = 'V'
public const GSVDQ: GSVDJob = 'Q'
public const GSVDUnit: GSVDJob = 'I'
public const GSVDNone: GSVDJob = 'N'

public type EVComp = UInt8

public const EVOrig: EVComp = 'V'
public const EVTridiag: EVComp = 'I'
public const EVCompNone: EVComp = 'N'

public type EVJob = UInt8

public const EVCompute: EVJob = 'V'
public const EVNone: EVJob = 'N'

public type LeftEVJob = UInt8

public const LeftEVCompute: LeftEVJob = 'V'
public const LeftEVNone: LeftEVJob = 'N'

public type RightEVJob = UInt8

public const RightEVCompute: RightEVJob = 'V'
public const RightEVNone: RightEVJob = 'N'

public type BalanceJob = UInt8

public const Permute: BalanceJob = 'P'
public const Scale: BalanceJob = 'S'
public const PermuteScale: BalanceJob = 'B'
public const BalanceNone: BalanceJob = 'N'

public type SchurJob = UInt8

public const EigenvaluesOnly: SchurJob = 'E'
public const EigenvaluesAndSchur: SchurJob = 'S'

public type SchurComp = UInt8

public const SchurOrig: SchurComp = 'V'
public const SchurHess: SchurComp = 'I'
public const SchurNone: SchurComp = 'N'

public type UpdateSchurComp = UInt8

public const UpdateSchur: UpdateSchurComp = 'V'
public const UpdateSchurNone: UpdateSchurComp = 'N'

public type EVSide = UInt8

public const EVRight: EVSide = 'R'
public const EVLeft: EVSide = 'L'
public const EVBoth: EVSide = 'B'

public type EVHowMany = UInt8

public const EVAll: EVHowMany = 'A'
public const EVAllMulQ: EVHowMany = 'B'
public const EVSelected: EVHowMany = 'S'

public type MaximizeNormXJob = UInt8

public const LocalLookAhead: MaximizeNormXJob = 0
public const NormalizedNullVector: MaximizeNormXJob = 2

public type OrthoComp = UInt8

public const OrthoNone: OrthoComp = 'N'
public const OrthoExplicit: OrthoComp = 'I'
public const OrthoPostmul: OrthoComp = 'V'
