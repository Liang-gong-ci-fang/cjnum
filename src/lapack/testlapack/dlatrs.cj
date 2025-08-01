package cjnum.lapack.testlapack

import std.math.*
import std.random.*
import std.collection.*
import cjnum.blas.*
import cjnum.blas.blas64.*
import cjnum.floats.*
import cjnum.lapack.*
import cjnum.util.*

public func dlatrsTest(lFloat64: LFloat64): Bool {
    let rnd = Random(1)
    for (uplo in [Upper, Lower]) {
        for (trans in [Trans, NoTrans]) {
            for (n in [0, 1, 2, 3, 4, 5, 6, 7, 10, 20, 50, 100]) {
                for (lda in [n, 2 * n + 1]) {
                    let lda_ = max(1, lda)
                    let imats = ArrayList<Int64>([7, 11, 12, 13, 14, 15, 16, 17, 18, 19])
                    if (n < 6) {
                        imats.add(19) 
                    }
                    for (imat in imats) {
                        _testDlatrs(lFloat64, imat, uplo, trans, n, lda_, rnd)
                    }
                }
            }
        }
    }
    return true
}

func _testDlatrs(lFloat64: LFloat64, imat: Int64, uplo: Uplo, trans: Transpose, n: Int64, lda: Int64, rnd: Random): Unit {
    const tol = 1e-14

    let a = nanSlice(n * lda)
    let b = nanSlice(n)
    let work = Array<Float64>(3 * n, repeat: 0.0)

    let diag = dlattr(imat, uplo, trans, n, a, lda, b, work, rnd)
    if (imat <= 10) {
        dlarnv(b, 3, rnd)
    }

    let cnorm = nanSlice(n)
    let x = b.clone()

    var scale = lFloat64.dlatrs(uplo, trans, diag, false, n, a, lda, x, cnorm)
    let prefix = "Case imat=${imat} (n=${n},lda=${lda},trans=${trans},uplo=${uplo},diag=${diag}"
    for (i in 0..cnorm.size) {
        if (cnorm[i].isNaN()) {
            throw Exception("${prefix}: cnorm[${i}] not computed (scale=${scale},normin=false)")
        }
    }
    var (resid, hasNaN) = dlatrsResidual(uplo, trans, diag, n, a, lda, scale, cnorm, x, b, work[..n])
    if (hasNaN) {
        throw Exception("${prefix}: unexpected NaN (scale=${scale},normin=false)")
    } else if (resid > tol) {
        throw Exception("${prefix}: residual ${resid} too large (scale=${scale},normin=false)")
    }

    copy(x, b)
    scale = lFloat64.dlatrs(uplo, trans, diag, true, n, a, lda, x, cnorm)
    (resid, hasNaN) = dlatrsResidual(uplo, trans, diag, n, a, lda, scale, cnorm, x, b, work[..n])
    if (hasNaN) {
        throw Exception("${prefix}: unexpected NaN (scale=${scale},normin=true)")
    } else if (resid > tol) {
        throw Exception("${prefix}: residual ${resid} too large (scale=${scale},normin=true)")
    }
}

func dlatrsResidual(uplo: Uplo, trans: Transpose, diag: Diag, n: Int64, a: Array<Float64>, lda: Int64, scale: Float64,
    cnorm: Array<Float64>, x: Array<Float64>, b: Array<Float64>, work: Array<Float64>): (Float64, Bool) {
    if (n == 0) {
        return (0.0, false)
    }

    var tnorm = 0.0
    if (diag == NonUnitDiag) {
        for (j in 0..n) {
            tnorm = nanMax(tnorm, abs(a[j * lda + j]) + cnorm[j])
        }
    } else {
        for (j in 0..n) {
            tnorm = nanMax(tnorm, 1.0 + cnorm[j])
        }
    }

    const eps = dlamchE
    const tiny = safmin

    let nFloat64 = nFloat64Implementation()
    copy(work, x)
    var ix = nFloat64.idamax(n, work, 1)
    var xnorm = nanMax(1.0, abs(work[ix]))
    let xscal = 1.0 / xnorm / Float64(n)
    nFloat64.dscal(n, xscal, work, 1)
    nFloat64.dtrmv(uplo, trans, diag, n, a, lda, work, 1)
    nFloat64.daxpy(n, -scale * xscal, b, 1, work, 1)
    for (v in work) {
        if (v.isNaN()) {
            return (1.0 / eps, true)
        }
    }
    ix = nFloat64.idamax(n, work, 1)
    var resid = abs(work[ix])
    ix = nFloat64.idamax(n, x, 1)
    xnorm = abs(x[ix])
    if (resid * tiny <= xnorm) {
        if (xnorm > 0.0) {
            resid /= xnorm
        }
    } else if (resid > 0.0) {
        resid = 1.0 / eps
    }
    if (resid * tiny <= tnorm) {
        if (tnorm > 0.0) {
            resid /= tnorm
        }
    } else if (resid > 0.0) {
        resid = 1.0 / eps
    }
    return (resid, false)
}

func copy(dst: Array<Float64>, src: Array<Float64>): Int64 {
    let minSize = min(src.size, dst.size)
    dst[0..minSize] = src[0..minSize]
    return minSize
}