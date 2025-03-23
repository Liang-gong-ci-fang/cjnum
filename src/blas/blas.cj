package cjnum.blas

import cjnum.complex.*

public type Flag = Int64

public const Identity: Flag = -2
public const Rescaling: Flag = -1
public const OffDiagonal: Flag = 0
public const Diagonal: Flag = 1

public class SrotmParams {
    public SrotmParams(public let flag: Flag, public let h: Array<Float32>) {}
}

public class DrotmParams {
    public DrotmParams(public let flag: Flag, public let h: Array<Float64>) {}
}

public type Transpose = UInt8

public const NoTrans: Transpose = 'N'
public const Trans: Transpose = 'T'
public const COnjTrans: Transpose = 'C'

public type Uplo = UInt8

public const Upper: Uplo = 'U'
public const Lower: Uplo = 'L'
public const All: Uplo = 'A'

public type Diag = UInt8

public const NonUnitDiag: Diag = 'N'
public const UnitDiag: Diag = 'U'

public type Side = UInt8

public const Left: Side = 'L'
public const Right: Side = 'R'

public interface NFloat32 <: Float32Level1 & Float32Level2 & Float32Level3 {}

public interface Float32Level1 {
    func sdsdot(n: Int64, alpha: Float32, x: Array<Float32>, incX: Int64, y: Array<Float32>, incY: Int64): Float32
    func dsdot(n: Int64, x: Array<Float32>, incX: Int64, y: Array<Float32>, incY: Int64): Float64
    func sdot(n: Int64, x: Array<Float32>, incX: Int64, y: Array<Float32>, incY: Int64): Float32
    func snrm2(n: Int64, x: Array<Float32>, incX: Int64): Float32
    func sasum(n: Int64, x: Array<Float32>, incX: Int64): Float32
    func isamax(n: Int64, x: Array<Float32>, incX: Int64): Int64
    func sswap(n: Int64, x: Array<Float32>, incX: Int64, y: Array<Float32>, incY: Int64): Unit
    func scopy(n: Int64, x: Array<Float32>, incX: Int64, y: Array<Float32>, incY: Int64): Unit
    func saxpy(n: Int64, alpha: Float32, x: Array<Float32>, incX: Int64, y: Array<Float32>, incY: Int64): Unit
    func srotg(a: Float32, b: Float32): (Float32, Float32, Float32, Float32)
    func srotmg(d1: Float32, d2: Float32, b1: Float32, b2: Float32): (SrotmParams, Float32, Float32, Float32)
    func srot(n: Int64, x: Array<Float32>, incX: Int64, y: Array<Float32>, incY: Int64, c: Float32, s: Float32): Unit
    func srotm(n: Int64, x: Array<Float32>, incX: Int64, y: Array<Float32>, incY: Int64, p: SrotmParams): Unit
    func sscal(n: Int64, alpha: Float32, x: Array<Float32>, incX: Int64): Unit
}

public interface Float32Level2 {
    func sgemv(tA: Transpose, m: Int64, n: Int64, alpha: Float32, a: Array<Float32>, lda: Int64, x: Array<Float32>,
        incX: Int64, beta: Float32, y: Array<Float32>, incY: Int64): Unit
    func sgbmv(tA: Transpose, m: Int64, n: Int64, kL: Int64, kU: Int64, alpha: Float32, a: Array<Float32>, lda: Int64,
        x: Array<Float32>, incX: Int64, beta: Float32, y: Array<Float32>, incY: Int64): Unit
    func strmv(ul: Uplo, tA: Transpose, d: Diag, n: Int64, a: Array<Float32>, lda: Int64, x: Array<Float32>, incX: Int64): Unit
    func stbmv(ul: Uplo, tA: Transpose, d: Diag, n: Int64, k: Int64, a: Array<Float32>, lda: Int64, x: Array<Float32>,
        incX: Int64): Unit
    func stpmv(ul: Uplo, tA: Transpose, d: Diag, n: Int64, ap: Array<Float32>, x: Array<Float32>, incX: Int64): Unit
    func strsv(ul: Uplo, tA: Transpose, d: Diag, n: Int64, a: Array<Float32>, lda: Int64, x: Array<Float32>, incX: Int64): Unit
    func stbsv(ul: Uplo, tA: Transpose, d: Diag, n: Int64, k: Int64, a: Array<Float32>, lda: Int64, x: Array<Float32>,
        incX: Int64): Unit
    func stpsv(ul: Uplo, tA: Transpose, d: Diag, n: Int64, ap: Array<Float32>, x: Array<Float32>, incX: Int64): Unit
    func ssymv(ul: Uplo, n: Int64, alpha: Float32, a: Array<Float32>, lda: Int64, x: Array<Float32>, incX: Int64,
        beta: Float32, y: Array<Float32>, incY: Int64): Unit
    func ssbmv(ul: Uplo, n: Int64, k: Int64, alpha: Float32, a: Array<Float32>, lda: Int64, x: Array<Float32>,
        incX: Int64, beta: Float32, y: Array<Float32>, incY: Int64): Unit
    func sspmv(ul: Uplo, n: Int64, alpha: Float32, ap: Array<Float32>, x: Array<Float32>, incX: Int64, beta: Float32,
        y: Array<Float32>, incY: Int64): Unit
    func sger(m: Int64, n: Int64, alpha: Float32, x: Array<Float32>, incX: Int64, y: Array<Float32>, incY: Int64,
        a: Array<Float32>, lda: Int64): Unit
    func ssyr(ul: Uplo, n: Int64, alpha: Float32, x: Array<Float32>, incX: Int64, a: Array<Float32>, lda: Int64): Unit
    func sspr(ul: Uplo, n: Int64, alpha: Float32, x: Array<Float32>, incX: Int64, ap: Array<Float32>): Unit
    func ssyr2(ul: Uplo, n: Int64, alpha: Float32, x: Array<Float32>, incX: Int64, y: Array<Float32>, incY: Int64,
        a: Array<Float32>, lda: Int64): Unit
    func sspr2(ul: Uplo, n: Int64, alpha: Float32, x: Array<Float32>, incX: Int64, y: Array<Float32>, incY: Int64,
        a: Array<Float32>, lda: Int64): Unit
}

public interface Float32Level3 {
    func sgemm(tA: Transpose, tB: Transpose, m: Int64, n: Int64, k: Int64, alpha: Float32, a: Array<Float32>,
        lda: Int64, b: Array<Float32>, ldb: Int64, beta: Float32, c: Array<Float32>, ldc: Int64): Unit
    func ssymm(s: Side, ul: Uplo, m: Int64, n: Int64, alpha: Float32, a: Array<Float32>, lda: Int64, b: Array<Float32>,
        ldb: Int64, beta: Float32, c: Array<Float32>, ldc: Int64): Unit
    func ssyrk(ul: Uplo, t: Transpose, n: Int64, k: Int64, alpha: Float32, a: Array<Float32>, lda: Int64, beta: Float32,
        c: Array<Float32>, ldc: Int64): Unit
    func ssyr2k(ul: Uplo, t: Transpose, n: Int64, k: Int64, alpha: Float32, a: Array<Float32>, lda: Int64,
        b: Array<Float32>, ldb: Int64, beta: Float32, c: Array<Float32>, ldc: Int64): Unit
    func strmm(s: Side, ul: Uplo, tA: Transpose, d: Diag, m: Int64, n: Int64, alpha: Float32, a: Array<Float32>,
        lda: Int64, b: Array<Float32>, ldb: Int64): Unit
    func strsm(s: Side, ul: Uplo, tA: Transpose, d: Diag, m: Int64, n: Int64, alpha: Float32, a: Array<Float32>,
        lda: Int64, b: Array<Float32>, ldb: Int64): Unit
}

public interface NFloat64 <: Float64Level1 & Float64Level2 & Float64Level3 {}

public interface Float64Level1 {
    func ddot(n: Int64, x: Array<Float64>, incX: Int64, y: Array<Float64>, incY: Int64): Float64
    func dnrm2(n: Int64, x: Array<Float64>, incX: Int64): Float64
    func dasum(n: Int64, x: Array<Float64>, incX: Int64): Float64
    func idamax(n: Int64, x: Array<Float64>, incX: Int64): Int64
    func dswap(n: Int64, x: Array<Float64>, incX: Int64, y: Array<Float64>, incY: Int64): Unit
    func dcopy(n: Int64, x: Array<Float64>, incX: Int64, y: Array<Float64>, incY: Int64): Unit
    func daxpy(n: Int64, alpha: Float64, x: Array<Float64>, incX: Int64, y: Array<Float64>, incY: Int64): Unit
    func drotg(a: Float64, b: Float64): (Float64, Float64, Float64)
    func drotmg(d1: Float64, d2: Float64, b1: Float64, b2: Float64): (DrotmParams, Float64, Float64, Float64)
    func drot(n: Int64, x: Array<Float64>, incX: Int64, y: Array<Float64>, incY: Int64, c: Float64, s: Float64): Unit
    func drotm(n: Int64, x: Array<Float64>, incX: Int64, y: Array<Float64>, incY: Int64, p: DrotmParams): Unit
    func dscal(n: Int64, alpha: Float64, x: Array<Float64>, incX: Int64): Unit
}

public interface Float64Level2 {
    func dgemv(tA: Transpose, m: Int64, n: Int64, alpha: Float64, a: Array<Float64>, lda: Int64, x: Array<Float64>,
        incX: Int64, beta: Float64, y: Array<Float64>, incY: Int64): Unit
    func dgbmv(tA: Transpose, m: Int64, n: Int64, kL: Int64, kU: Int64, alpha: Float64, a: Array<Float64>, lda: Int64,
        x: Array<Float64>, incX: Int64, beta: Float64, y: Array<Float64>, incY: Int64): Unit
    func dtrmv(ul: Uplo, tA: Transpose, d: Diag, n: Int64, a: Array<Float64>, lda: Int64, x: Array<Float64>, incX: Int64): Unit
    func dtbmv(ul: Uplo, tA: Transpose, d: Diag, n: Int64, k: Int64, a: Array<Float64>, lda: Int64, x: Array<Float64>,
        incX: Int64): Unit
    func dtpmv(ul: Uplo, tA: Transpose, d: Diag, n: Int64, ap: Array<Float64>, x: Array<Float64>, incX: Int64): Unit
    func dtrsv(ul: Uplo, tA: Transpose, d: Diag, n: Int64, a: Array<Float64>, lda: Int64, x: Array<Float64>, incX: Int64): Unit
    func dtbsv(ul: Uplo, tA: Transpose, d: Diag, n: Int64, k: Int64, a: Array<Float64>, lda: Int64, x: Array<Float64>,
        incX: Int64): Unit
    func dtpsv(ul: Uplo, tA: Transpose, d: Diag, n: Int64, ap: Array<Float64>, x: Array<Float64>, incX: Int64): Unit
    func dsymv(ul: Uplo, n: Int64, alpha: Float64, a: Array<Float64>, lda: Int64, x: Array<Float64>, incX: Int64,
        beta: Float64, y: Array<Float64>, incY: Int64): Unit
    func dsbmv(ul: Uplo, n: Int64, k: Int64, alpha: Float64, a: Array<Float64>, lda: Int64, x: Array<Float64>,
        incX: Int64, beta: Float64, y: Array<Float64>, incY: Int64): Unit
    func dspmv(ul: Uplo, n: Int64, alpha: Float64, ap: Array<Float64>, x: Array<Float64>, incX: Int64, beta: Float64,
        y: Array<Float64>, incY: Int64): Unit
    func dger(m: Int64, n: Int64, alpha: Float64, x: Array<Float64>, incX: Int64, y: Array<Float64>, incY: Int64,
        a: Array<Float64>, lda: Int64): Unit
    func dsyr(ul: Uplo, n: Int64, alpha: Float64, x: Array<Float64>, incX: Int64, a: Array<Float64>, lda: Int64): Unit
    func dspr(ul: Uplo, n: Int64, alpha: Float64, x: Array<Float64>, incX: Int64, ap: Array<Float64>): Unit
    func dsyr2(ul: Uplo, n: Int64, alpha: Float64, x: Array<Float64>, incX: Int64, y: Array<Float64>, incY: Int64,
        a: Array<Float64>, lda: Int64): Unit
    func dspr2(ul: Uplo, n: Int64, alpha: Float64, x: Array<Float64>, incX: Int64, y: Array<Float64>, incY: Int64,
        a: Array<Float64>): Unit
}

public interface Float64Level3 {
    func dgemm(tA: Transpose, tB: Transpose, m: Int64, n: Int64, k: Int64, alpha: Float64, a: Array<Float64>,
        lda: Int64, b: Array<Float64>, ldb: Int64, beta: Float64, c: Array<Float64>, ldc: Int64): Unit
    func dsymm(s: Side, ul: Uplo, m: Int64, n: Int64, alpha: Float64, a: Array<Float64>, lda: Int64, b: Array<Float64>,
        ldb: Int64, beta: Float64, c: Array<Float64>, ldc: Int64): Unit
    func dsyrk(ul: Uplo, t: Transpose, n: Int64, k: Int64, alpha: Float64, a: Array<Float64>, lda: Int64, beta: Float64,
        c: Array<Float64>, ldc: Int64): Unit
    func dsyr2k(ul: Uplo, t: Transpose, n: Int64, k: Int64, alpha: Float64, a: Array<Float64>, lda: Int64,
        b: Array<Float64>, ldb: Int64, beta: Float64, c: Array<Float64>, ldc: Int64): Unit
    func dtrmm(s: Side, ul: Uplo, tA: Transpose, d: Diag, m: Int64, n: Int64, alpha: Float64, a: Array<Float64>,
        lda: Int64, b: Array<Float64>, ldb: Int64): Unit
    func dtrsm(s: Side, ul: Uplo, tA: Transpose, d: Diag, m: Int64, n: Int64, alpha: Float64, a: Array<Float64>,
        lda: Int64, b: Array<Float64>, ldb: Int64): Unit
}

public interface Complex64 <: Complex64Level1 & Complex64Level2 & Complex64Level3 {}

public interface Complex64Level1 {
    func cdotu(n: Int64, x: Array<Complex64>, incX: Int64, y: Array<Complex64>, incY: Int64): Complex64
    func cdotc(n: Int64, x: Array<Complex64>, incX: Int64, y: Array<Complex64>, incY: Int64): Complex64
    func scnrm2(n: Int64, x: Array<Complex64>, incX: Int64): Float32
    func scasum(n: Int64, x: Array<Complex64>, incX: Int64): Float32
    func icamax(n: Int64, x: Array<Complex64>, incX: Int64): Int64
    func cswap(n: Int64, x: Array<Complex64>, incX: Int64, y: Array<Complex64>, incY: Int64): Unit
    func ccopy(n: Int64, x: Array<Complex64>, incX: Int64, y: Array<Complex64>, incY: Int64): Unit
    func caxpy(n: Int64, alpha: Complex64, x: Array<Complex64>, incX: Int64, y: Array<Complex64>, incY: Int64): Unit
    func cscal(n: Int64, alpha: Complex64, x: Array<Complex64>, incX: Int64): Unit
    func csscal(n: Int64, alpha: Float32, x: Array<Complex64>, incX: Int64): Unit
}

public interface Complex64Level2 {
    func cgemv(tA: Transpose, m: Int64, n: Int64, alpha: Complex64, a: Array<Complex64>, lda: Int64,
        x: Array<Complex64>, incX: Int64, beta: Complex64, y: Array<Complex64>, incY: Int64): Unit
    func cgbmv(tA: Transpose, m: Int64, n: Int64, kL: Int64, alpha: Complex64, a: Array<Complex64>, lda: Int64,
        x: Array<Complex64>, incX: Int64, beta: Complex64, y: Array<Complex64>, incY: Int64): Unit
    func ctrmv(ul: Uplo, tA: Transpose, d: Diag, n: Int64, a: Array<Complex64>, lda: Int64, x: Array<Complex64>,
        incX: Int64): Unit
    func ctbmv(ul: Uplo, tA: Transpose, d: Diag, n: Int64, k: Int64, a: Array<Complex64>, lda: Int64,
        x: Array<Complex64>, incX: Int64): Unit
    func ctpmv(ul: Uplo, tA: Transpose, d: Diag, n: Int64, ap: Array<Complex64>, x: Array<Complex64>, incX: Int64): Unit
    func ctrsv(ul: Uplo, tA: Transpose, d: Diag, n: Int64, a: Array<Complex64>, lda: Int64, x: Array<Complex64>,
        incX: Int64): Unit
    func ctbsv(ul: Uplo, tA: Transpose, d: Diag, n: Int64, k: Int64, a: Array<Complex64>, lda: Int64,
        x: Array<Complex64>, incX: Int64): Unit
    func ctpsv(ul: Uplo, tA: Transpose, d: Diag, n: Int64, ap: Array<Complex64>, x: Array<Complex64>, incX: Int64): Unit
    func chemv(ul: Uplo, n: Int64, alpha: Complex64, a: Array<Complex64>, lda: Int64, x: Array<Complex64>, incX: Int64,
        beta: Complex64, y: Array<Complex64>, incY: Int64): Unit
    func chbmv(ul: Uplo, n: Int64, k: Int64, alpha: Complex64, a: Array<Complex64>, lda: Int64, x: Array<Complex64>,
        incX: Int64, beta: Complex64, y: Array<Complex64>, incY: Int64): Unit
    func chpmv(ul: Uplo, n: Int64, alpha: Complex64, ap: Array<Complex64>, x: Array<Complex64>, incX: Int64,
        beta: Complex64, y: Array<Complex64>, incY: Int64): Unit
    func cgeru(m: Int64, n: Int64, alpha: Complex64, x: Array<Complex64>, incX: Int64, y: Array<Complex64>, incY: Int64,
        a: Array<Complex64>, lda: Int64): Unit
    func cgerc(m: Int64, n: Int64, alpha: Complex64, x: Array<Complex64>, incX: Int64, y: Array<Complex64>, incY: Int64,
        a: Array<Complex64>, lda: Int64): Unit
    func cher(ul: Uplo, n: Int64, alpha: Float32, x: Array<Complex64>, incX: Int64, a: Array<Complex64>, lda: Int64): Unit
    func chpr(ul: Uplo, n: Int64, alpha: Float32, x: Array<Complex64>, incX: Int64, a: Array<Complex64>): Unit
    func cher2(ul: Uplo, n: Int64, alpha: Complex64, x: Array<Complex64>, incX: Int64, y: Array<Complex64>, incY: Int64,
        a: Array<Complex64>, lda: Int64): Unit
    func chpr2(ul: Uplo, n: Int64, alpha: Complex64, x: Array<Complex64>, incX: Int64, y: Array<Complex64>,
        a: Array<Complex64>, lda: Int64): Unit
}

public interface Complex64Level3 {
    func cgemm(tA: Transpose, tB: Transpose, m: Int64, n: Int64, k: Int64, alpha: Complex64, a: Array<Complex64>,
        lda: Int64, b: Array<Complex64>, ldb: Int64, beta: Complex64, c: Array<Complex64>, ldc: Int64): Unit
    func csymm(side: Side, ul: Uplo, m: Int64, n: Int64, alpha: Complex64, a: Array<Complex64>, lda: Int64,
        b: Array<Complex64>, ldb: Int64, beta: Complex64, c: Array<Complex64>, ldc: Int64): Unit
    func csyrk(ul: Uplo, t: Transpose, n: Int64, k: Int64, alpha: Complex64, a: Array<Complex64>, lda: Int64,
        beta: Complex64, c: Array<Complex64>, ldc: Int64): Unit
    func csyr2k(ul: Uplo, t: Transpose, n: Int64, k: Int64, alpha: Complex64, a: Array<Complex64>, lda: Int64,
        b: Array<Complex64>, ldb: Int64, beta: Complex64, c: Array<Complex64>, ldc: Int64): Unit
    func ctrmm(side: Side, ul: Uplo, tA: Transpose, d: Diag, m: Int64, n: Int64, alpha: Complex64, a: Array<Complex64>,
        lda: Int64, b: Array<Complex64>, ldb: Int64): Unit
    func ctrsm(side: Side, ul: Uplo, tA: Transpose, d: Diag, m: Int64, n: Int64, alpha: Complex64, a: Array<Complex64>,
        lda: Int64, b: Array<Complex64>, ldb: Int64): Unit
    func chemm(side: Side, ul: Uplo, m: Int64, n: Int64, alpha: Complex64, a: Array<Complex64>, lda: Int64,
        b: Array<Complex64>, ldb: Int64, beta: Complex64, c: Array<Complex64>, ldc: Int64): Unit
    func cherk(ul: Uplo, t: Transpose, n: Int64, k: Int64, alpha: Float32, a: Array<Complex64>, lda: Int64,
        beta: Float32, c: Array<Complex64>, ldc: Int64): Unit
    func cher2k(ul: Uplo, t: Transpose, n: Int64, k: Int64, alpha: Complex64, a: Array<Complex64>, lda: Int64,
        b: Array<Complex64>, ldb: Int64, beta: Float32, c: Array<Complex64>, ldc: Int64): Unit
}

public interface NComplex128 <: Complex128Level1 & Complex128Level2 & Complex128Level3 {}

public interface Complex128Level1 {
    func zdotu(n: Int64, x: Array<Complex128>, incX: Int64, y: Array<Complex128>, incY: Int64): Complex128
    func zdotc(n: Int64, x: Array<Complex128>, incX: Int64, y: Array<Complex128>, incY: Int64): Complex128
    func dzrm2(n: Int64, x: Array<Complex128>, incX: Int64): Float64
    func dzasum(n: Int64, x: Array<Complex128>, incX: Int64): Float64
    func izamax(n: Int64, x: Array<Complex128>, incX: Int64): Int64
    func zswap(n: Int64, x: Array<Complex128>, incX: Int64, y: Array<Complex128>, incY: Int64): Unit
    func zcopy(n: Int64, x: Array<Complex128>, incX: Int64, y: Array<Complex128>, incY: Int64): Unit
    func zaxpy(n: Int64, alpha: Complex128, x: Array<Complex128>, incX: Int64, y: Array<Complex128>, incY: Int64): Unit
    func zscal(n: Int64, alpha: Complex128, x: Array<Complex128>, incX: Int64): Unit
    func zdscal(n: Int64, alpha: Float64, x: Array<Complex128>, incX: Int64): Unit
}

public interface Complex128Level2 {
    func zgemv(tA: Transpose, m: Int64, n: Int64, alpha: Complex128, a: Array<Complex128>, lda: Int64,
        x: Array<Complex128>, incX: Int64, beta: Complex128, y: Array<Complex128>, incY: Int64): Unit
    func zgbmv(tA: Transpose, m: Int64, n: Int64, kL: Int64, kU: Int64, alpha: Complex128, a: Array<Complex128>,
        lda: Int64, x: Array<Complex128>, incX: Int64, beta: Complex128, y: Array<Complex128>, incY: Int64): Unit
    func ztrmv(ul: Uplo, tA: Transpose, d: Diag, n: Int64, a: Array<Complex128>, lda: Int64, x: Array<Complex128>,
        incX: Int64): Unit
    func ztbmv(ul: Uplo, tA: Transpose, d: Diag, n: Int64, k: Int64, a: Array<Complex128>, lda: Int64,
        x: Array<Complex128>, incX: Int64): Unit
    func ztpmv(ul: Uplo, tA: Transpose, d: Diag, n: Int64, ap: Array<Complex128>, x: Array<Complex128>, incX: Int64): Unit
    func ztrsv(ul: Uplo, tA: Transpose, d: Diag, n: Int64, a: Array<Complex128>, lda: Int64, x: Array<Complex128>,
        incX: Int64): Unit
    func ztbsv(ul: Uplo, tA: Transpose, d: Diag, n: Int64, k: Int64, a: Array<Complex128>, lda: Int64,
        x: Array<Complex128>, incX: Int64): Unit
    func ztpsv(ul: Uplo, tA: Transpose, d: Diag, n: Int64, ap: Array<Complex128>, x: Array<Complex128>, incX: Int64): Unit
    func zhemv(ul: Uplo, n: Int64, alpha: Complex128, a: Array<Complex128>, lda: Int64, x: Array<Complex128>,
        incX: Int64, beta: Complex128, y: Array<Complex128>, incY: Int64): Unit
    func zhbmv(ul: Uplo, n: Int64, k: Int64, alpha: Complex128, a: Array<Complex128>, lda: Int64, x: Array<Complex128>,
        incX: Int64, beta: Complex128, y: Array<Complex128>, incY: Int64): Unit
    func zhpmv(ul: Uplo, n: Int64, alpha: Complex128, ap: Array<Complex128>, x: Array<Complex128>, incX: Int64,
        beta: Complex128, y: Array<Complex128>, incY: Int64): Unit
    func zgeru(m: Int64, n: Int64, alpha: Complex128, x: Array<Complex128>, incX: Int64, y: Array<Complex128>,
        incY: Int64, a: Array<Complex128>, lda: Int64): Unit
    func zgerc(m: Int64, n: Int64, alpha: Complex128, x: Array<Complex128>, incX: Int64, y: Array<Complex128>,
        incY: Int64, a: Array<Complex128>, lda: Int64): Unit
    func zher(ul: Uplo, n: Int64, alpha: Float64, x: Array<Complex128>, incX: Int64, a: Array<Complex128>, lda: Int64): Unit
    func zhpr(ul: Uplo, n: Int64, alpha: Float64, x: Array<Complex128>, incX: Int64, a: Array<Complex128>): Unit
    func zher2(ul: Uplo, n: Int64, alpha: Complex128, x: Array<Complex128>, incX: Int64, y: Array<Complex128>,
        incY: Int64, a: Array<Complex128>, lda: Int64): Unit
    func zhpr2(ul: Uplo, n: Int64, alpha: Complex128, x: Array<Complex128>, incX: Int64, y: Array<Complex128>,
        incY: Int64, ap: Array<Complex128>): Unit
}

public interface Complex128Level3 {
    func zgemm(tA: Transpose, tB: Transpose, m: Int64, n: Int64, k: Int64, alpha: Complex128, a: Array<Complex128>,
        lda: Int64, b: Array<Complex128>, ldb: Int64, beta: Complex128, c: Array<Complex128>, ldc: Int64): Unit
    func zsymm(side: Side, ul: Uplo, m: Int64, n: Int64, alpha: Complex128, a: Array<Complex128>, lda: Int64,
        b: Array<Complex128>, ldb: Int64, beta: Complex128, c: Array<Complex128>, ldc: Int64): Unit
    func zsyrk(ul: Uplo, t: Transpose, n: Int64, k: Int64, alpha: Complex128, a: Array<Complex128>, lda: Int64,
        beta: Complex128, c: Array<Complex128>, ldc: Int64): Unit
    func zsyr2k(ul: Uplo, t: Transpose, n: Int64, k: Int64, alpha: Complex128, a: Array<Complex128>, lda: Int64,
        b: Array<Complex128>, ldb: Int64, beta: Complex128, c: Array<Complex128>, ldc: Int64): Unit
    func ztrmm(side: Side, ul: Uplo, tA: Transpose, d: Diag, m: Int64, n: Int64, alpha: Complex128,
        a: Array<Complex128>, lda: Int64, b: Array<Complex128>, ldb: Int64): Unit
    func ztrsm(side: Side, ul: Uplo, tA: Transpose, d: Diag, m: Int64, n: Int64, alpha: Complex128,
        a: Array<Complex128>, lda: Int64, b: Array<Complex128>, ldb: Int64): Unit
    func zhemm(side: Side, ul: Uplo, m: Int64, n: Int64, alpha: Complex128, a: Array<Complex128>, lda: Int64,
        b: Array<Complex128>, ldb: Int64, beta: Complex128, c: Array<Complex128>, ldc: Int64): Unit
    func zherk(ul: Uplo, t: Transpose, n: Int64, k: Int64, alpha: Float64, a: Array<Complex128>, lda: Int64,
        beta: Float64, c: Array<Complex128>, ldc: Int64): Unit
    func zher2k(ul: Uplo, t: Transpose, n: Int64, k: Int64, alpha: Complex128, a: Array<Complex128>, lda: Int64,
        b: Array<Complex128>, ldb: Int64, beta: Float64, c: Array<Complex128>, ldc: Int64): Unit
}
