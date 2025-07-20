package cjnum.blas.blas64

import cjnum.blas.cjnum.*
import cjnum.blas.*
import cjnum.*
import cjnum.blas.*
import cjnum.blas.cjnum.*

var blas64: NFloat64 = NFloat64Implementation()

public func nFloat64use(b: NFloat64): Unit {
    blas64 = b
}

public func nFloat64Implementation(): NFloat64 {
    return blas64
}

public func flatten(a: Array<Array<Float64>>): Array<Float64> {
    if (a.size == 0) {
        return []
    }
    let m = a.size
    let n = a[0].size
    let s = Array<Float64>(m * n, repeat: 0.0)
    for (i in 0..m) {
        for (j in 0..n) {
            s[i * n + j] = a[i][j]
        }
    }
    return s
}

private const negInc    = "blas64: negative vector increment"
private const badLength = "blas64: vector length mismatch"

public func dot(x: Vector, y: Vector): Float64 {
    if (x.n != y.n) {
        throw Exception(badLength)
    }
    return blas64.ddot(x.n, x.data, x.inc, y.data, y.inc)
}

public func nrm2(x: Vector): Float64 {
    if (x.inc < 0) {
        throw Exception(negInc)
    }
    return blas64.dnrm2(x.n, x.data, x.inc)
}

public func asum(x: Vector): Float64 {
    if (x.inc < 0) {
        throw Exception(negInc)
    }
    return blas64.dasum(x.n, x.data, x.inc)
}

public func iamax(x: Vector): Int64 {
    if (x.inc < 0) {
        throw Exception(negInc)
    }
    return blas64.idamax(x.n, x.data, x.inc)
}

public func swap(x: Vector, y: Vector): Unit {
    if (x.n != y.n) {
        throw Exception(badLength)
    } 
    blas64.dswap(x.n, x.data, x.inc, y.data, y.inc)
}

public func copy(x: Vector, y: Vector): Unit {
    if (x.n != y.n) {
        throw Exception(badLength)
    }
    blas64.dcopy(x.n, x.data, x.inc, y.data, y.inc)
}

public func axpy(alpha: Float64, x: Vector, y: Vector): Unit {
    if (x.n != y.n) {
        throw Exception(badLength)
    }
    blas64.daxpy(x.n, alpha, x.data, x.inc, y.data, y.inc)
}

public func Rotq(a: Float64, b: Float64): (Float64, Float64, Float64, Float64) {
    return blas64.drotg(a, b)
}

public func rotmq(d1: Float64, d2: Float64, b1: Float64, b2: Float64): (DrotmParams, Float64, Float64, Float64) {
    return blas64.drotmg(d1, d2, b1, b2)
}

public func rot(x: Vector, y: Vector, c: Float64, s: Float64): Unit {
    if (x.n != y.n) {
        throw Exception(badLength)
    }
    blas64.drot(x.n, x.data, x.inc, y.data, y.inc, c, s)
}

public func rotm(x: Vector, y: Vector, p: DrotmParams): Unit {
    if (x.n != y.n) {
        throw Exception(badLength)
    }
    blas64.drotm(x.n, x.data, x.inc, y.data, y.inc, p)
}

public func scal(alpha: Float64, x: Vector): Unit {
    if (x.inc < 0) {
        throw Exception(negInc)
    }
    blas64.dscal(x.n, alpha, x.data, x.inc)
}

public func gemv(t: Transpose, alpha: Float64, a: General, x: Vector, beta: Float64, y: Vector): Unit {
    blas64.dgemv(t, a.rows, a.cols, alpha, a.data, a.stride, x.data, x.inc, beta, y.data, y.inc)
}

public func gbmv(t: Transpose, alpha: Float64, a: Band, x: Vector, beta: Float64, y: Vector): Unit {
    blas64.dgbmv(t, a.rows, a.cols, a.kl, a.ku, alpha, a.data, a.stride, x.data, x.inc, beta, y.data, y.inc)
}

public func trmv(t: Transpose, a: Triangular, x: Vector): Unit {
    blas64.dtrmv(a.uplo, t, a.diag, a.n, a.data, a.stride, x.data, x.inc)
}

public func tbmv(t: Transpose, a: TriangularBand, x: Vector): Unit {
    blas64.dtbmv(a.uplo, t, a.diag, a.n, a.k, a.data, a.stride, x.data, x.inc)
}

public func tpmv(t: Transpose, a: TriangularPacked, x: Vector): Unit {
    blas64.dtpmv(a.uplo, t, a.diag, a.n, a.data, x.data, x.inc)
}

public func trsv(t: Transpose, a: Triangular, x: Vector): Unit {
    blas64.dtrsv(a.uplo, t, a.diag, a.n, a.data, a.stride, x.data, x.inc)
}

public func tbsv(t: Transpose, a: TriangularBand, x: Vector): Unit {
    blas64.dtbsv(a.uplo, t, a.diag, a.n, a.k, a.data, a.stride, x.data, x.inc)
}

public func tpsv(t: Transpose, a: TriangularPacked, x: Vector): Unit {
    blas64.dtpsv(a.uplo, t, a.diag, a.n, a.data, x.data, x.inc)
}

public func symv(alpha: Float64, a: Symmetric, x: Vector, beta: Float64, y: Vector): Unit {
    blas64.dsymv(a.uplo, a.n, alpha, a.data, a.stride, x.data, x.inc, beta, y.data, y.inc)
}

public func sbmv(alpha: Float64, a: SymmetricBand, x: Vector, beta: Float64, y: Vector): Unit {
    blas64.dsbmv(a.uplo, a.n, a.k, alpha, a.data, a.stride, x.data, x.inc, beta, y.data, y.inc)
}

public func spmv(alpha: Float64, a: SymmetricPacked, x: Vector, beta: Float64, y: Vector): Unit {
    blas64.dspmv(a.uplo, a.n, alpha, a.data, x.data, x.inc, beta, y.data, y.inc)
}

public func ger(alpha: Float64, x: Vector, y: Vector, a: General): Unit {
    blas64.dger(a.rows, a.cols, alpha, x.data, x.inc, y.data, y.inc, a.data, a.stride)
}

public func syr(alpha: Float64, x: Vector, a: Symmetric): Unit {
    blas64.dsyr(a.uplo, a.n, alpha, x.data, x.inc, a.data, a.stride)
}

public func spr(alpha: Float64, x: Vector, a: SymmetricPacked): Unit {
    blas64.dspr(a.uplo, a.n, alpha, x.data, x.inc, a.data)
}

public func syr2(alpha: Float64, x: Vector, y: Vector, a: Symmetric): Unit {
    blas64.dsyr2(a.uplo, a.n, alpha, x.data, x.inc, y.data, y.inc, a.data, a.stride)
}

public func spr2(alpha: Float64, x: Vector, y: Vector, a: SymmetricPacked): Unit {
    blas64.dspr2(a.uplo, a.n, alpha, x.data, x.inc, y.data, y.inc, a.data)
}

public func gemm(tA: Transpose, tB: Transpose, alpha: Float64, a: General, b: General, beta: Float64, c: General): Unit {
    let (m, k) = if (tA == NoTrans) {
        (a.rows, a.cols)
    } else {
        (a.cols, a.rows)
    }
    let n = if (tB == NoTrans) {
        b.cols
    } else {
        b.rows
    }
    blas64.dgemm(tA, tB, m, n, k, alpha, a.data, a.stride, b.data, b.stride, beta, c.data, c.stride)
}

public func symm(s: Side, alpha: Float64, a: Symmetric, b: General, beta: Float64, c: General): Unit {
    let (m, n) = if (s == Left) {
        (a.n, b.cols)
    } else {
        (b.rows, a.n)
    }
    blas64.dsymm(s, a.uplo, m, n, alpha, a.data, a.stride, b.data, b.stride, beta, c.data, c.stride)
}

public func syrk(t: Transpose, alpha: Float64, a: General, beta: Float64, c: Symmetric): Unit {
    let (n, k) = if (t == NoTrans) {
        (a.rows, a.cols)
    } else {
        (a.cols, a.rows)
    }
    blas64.dsyrk(c.uplo, t, n, k, alpha, a.data, a.stride, beta, c.data, c.stride)
}

public func syr2k(t: Transpose, alpha: Float64, a: General, b: General, beta: Float64, c: Symmetric): Unit {
    let (n, k) = if (t == NoTrans) {
        (a.rows, a.cols)
    } else {
        (a.cols, a.rows)
    }
    blas64.dsyr2k(c.uplo, t, n, k, alpha, a.data, a.stride, b.data, b.stride, beta, c.data, c.stride)
}

public func trmm(s: Side, tA: Transpose, alpha: Float64, a: Triangular, b: General): Unit {
    blas64.dtrmm(s, a.uplo, tA, a.diag, b.rows, b.cols, alpha, a.data, a.stride, b.data, b.stride)
}

public func trsm(s: Side, tA: Transpose, alpha: Float64, a: Triangular, b: General): Unit {
    blas64.dtrsm(s, a.uplo, tA, a.diag, b.rows, b.cols, alpha, a.data, a.stride, b.data, b.stride)
}