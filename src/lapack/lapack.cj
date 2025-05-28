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
    func dlasy2(tranl: Bool, tranr: Bool, isgn: Int, n1: Int, n2: Int, tl: Array<Float64>, ldtl: Int,
        tr: Array<Float64>, ldtr: Int, b: Array<Float64>, ldb: Int, x: Array<Float64>, ldx: Int): (Float64, Float64, 
        Bool)
    func dlascl(kind: MatrixType, kl: Int, ku: Int, cfrom: Float64, cto: Float64, m: Int, n: Int, a: Array<Float64>,
        lda: Int): Unit
    func dlaexc(wantq: Bool, n: Int, t: Array<Float64>, ldt: Int, q: Array<Float64>, ldq: Int, j1: Int, n1: Int,
        n2: Int, work: Array<Float64>): Bool
    func dlasv2(f: Float64, g: Float64, h: Float64): (ssmin: Float64, ssmax: Float64, snr: Float64, csr: Float64, 
        snl: Float64, csl: Float64)
    func dlauu2(uplo: Uplo, n: Int, a: Array<Float64>, lda: Int): Unit
    func dpotf2(uplo: Uplo, n: Int, a: Array<Float64>, lda: Int): Bool
    func dpbtf2(uplo: Uplo, n: Int, kd: Int, ab: Array<Float64>, ldab: Int): Bool
    func dgetc2(n: Int, a: Array<Float64>, lda: Int, ipiv: Array<Int>, jpiv: Array<Int>): Int
    func dgesc2(n: Int, a: Array<Float64>, lda: Int, rhs: Array<Float64>, ipiv: Array<Int>, jpiv: Array<Int>): Float64
    func dgghrd(compq: OrthoComp, compz: OrthoComp, n: Int, ilo: Int, ihi: Int, a: Array<Float64>, lda: Int,
        b: Array<Float64>, ldb: Int, q: Array<Float64>, ldq: Int, z: Array<Float64>, ldz: Int): Unit
}

type Direct = UInt8

public const Forward: Direct = 'F'
public const Backward: Direct = 'B'

type Sort = UInt8

public const SortIncreasing: Sort = 'I'
public const SortDecreasing: Sort = 'D'

type StoreV = UInt8

public const ColumnWise: StoreV = 'C'
public const RowWise: StoreV = 'R'

type MatrixNorm = UInt8

public const MaxAbs: MatrixNorm = 'M'
public const MaxColumnSum: MatrixNorm = 'O'
public const MaxRowSum: MatrixNorm = 'I'
public const Frobenius: MatrixNorm = 'F'

type MatrixType = UInt8

public const General: MatrixType = 'G'
public const UpperTri: MatrixType = 'U'
public const LowerTri: MatrixType = 'L'

type Pivot = UInt8

public const Variable: Pivot = 'V'
public const Top: Pivot = 'T'
public const Bottom: Pivot = 'B'

type ApplyOrtho = UInt8

public const ApplyP: ApplyOrtho = 'P'
public const ApplyQ: ApplyOrtho = 'Q'

type GenOrtho = UInt8

public const GeneratePT: GenOrtho = 'P'
public const GenerateQ: GenOrtho = 'Q'

type SVDJob = UInt8

public const SVDAll: SVDJob = 'A'
public const SVDStore: SVDJob = 'S'
public const SVDOverwrite: SVDJob = 'O'
public const SVDNone: SVDJob = 'N'

type GSVDJob = UInt8

public const GSVDU: GSVDJob = 'U'
public const GSVDV: GSVDJob = 'V'
public const GSVDQ: GSVDJob = 'Q'
public const GSVDUnit: GSVDJob = 'I'
public const GSVDNone: GSVDJob = 'N'

type EVComp = UInt8

public const EVOrig: EVComp = 'V'
public const EVTridiag: EVComp = 'I'
public const EVCompNone: EVComp = 'N'

type EVJob = UInt8

public const EVCompute: EVJob = 'V'
public const EVNone: EVJob = 'N'

type LeftEVJob = UInt8

public const LeftEVCompute: LeftEVJob = 'V'
public const LeftEVNone: LeftEVJob = 'N'

type RightEVJob = UInt8

public const RightEVCompute: RightEVJob = 'V'
public const RightEVNone: RightEVJob = 'N'

type BalanceJob = UInt8

public const Permute: BalanceJob = 'P'
public const Scale: BalanceJob = 'S'
public const PermuteScale: BalanceJob = 'B'
public const BalanceNone: BalanceJob = 'N'

type SchurJob = UInt8

public const EigenvaluesOnly: SchurJob = 'E'
public const EigenvaluesAndSchur: SchurJob = 'S'

type SchurComp = UInt8

public const SchurOrig: SchurComp = 'V'
public const SchurHess: SchurComp = 'I'
public const SchurNone: SchurComp = 'N'

type UpdateSchurComp = UInt8

public const UpdateSchur: UpdateSchurComp = 'V'
public const UpdateSchurNone: UpdateSchurComp = 'N'

type EVSide = UInt8

public const EVRight: EVSide = 'R'
public const EVLeft: EVSide = 'L'
public const EVBoth: EVSide = 'B'

type EVHowMany = UInt8

public const EVAll: EVHowMany = 'A'
public const EVAllMulQ: EVHowMany = 'B'
public const EVSelected: EVHowMany = 'S'

type MaximizeNormXJob = UInt8

public const LocalLookAhead: MaximizeNormXJob = 0
public const NormalizedNullVector: MaximizeNormXJob = 2

type OrthoComp = UInt8

public const OrthoNone: OrthoComp = 'N'
public const OrthoExplicit: OrthoComp = 'I'
public const OrthoPostmul: OrthoComp = 'V'
