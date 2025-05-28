package cjnum.lapack.testlapack

import cjnum.blas.blas64.*
import cjnum.blas.*
import cjnum.floats.*
import cjnum.lapack.*
import std.math.*
import std.random.*

public func dgesc2Test(lFloat64: LFloat64): Bool {
    let rnd = Random(1)
    for (n in [0, 1, 2, 3, 4, 5, 10, 20, 50]) {
        for (lda in [n, n + 3]) {
            if (!testDgesc2(lFloat64, rnd, n, lda, false)) {
                return false
            }
            if (!testDgesc2(lFloat64, rnd, n, lda, true)) {
                return false
            }
        }
    }
    return true
}

func testDgesc2(lFloat64: LFloat64, rnd: Random, n: Int, lda: Int, big: Bool): Bool {
    const tol = 1e-14

    let name = "n=${n},lda=${lda},big=${big}"

    let a = randomGeneral(n, n, max(1, lda), rnd)

    let b = randomGeneral(n, 1, 1, rnd)
    if (big) {
        for (i in 0..n) {
            b.data[i] *= bignum
        }
    }

    let lu = cloneGeneral(a)
    let ipiv = Array<Int>(n, item: 0)
    let jpiv = Array<Int>(n, item: 0)
    lFloat64.dgetc2(n, lu.data, lu.stride, ipiv, jpiv)

    let luCopy = cloneGeneral(lu)
    let ipivCopy = Array<Int>(ipiv.size, item: 0)
    copy(ipivCopy, ipiv)
    let jpivCopy = Array<Int>(jpiv.size, item: 0)
    copy(jpivCopy, jpiv)

    let x = cloneGeneral(b)
    let scale = lFloat64.dgesc2(n, lu.data, lu.stride, x.data, ipiv, jpiv)

    if (n == 0) {
        return true
    }

    if (!fSame(lu.data, luCopy.data)) {
        throw Exception("${name}: unexpected modification in lu")
    }
    if (!intsEqual(ipiv, ipivCopy)) {
        throw Exception("${name}: unexpected modification in ipiv")
    }
    if (!intsEqual(jpiv, jpivCopy)) {
        throw Exception("${name}: unexpected modification in jpiv")
    }

    if (scale <= 0.0 || 1.0 < scale) {
        throw Exception("${name}: scale ${scale} out of bounds (0,1]")
    }
    if (!big && scale != 1.0) {
        throw Exception("${name}: unexpected scaling, scale=${scale}")
    }

    let diff = b
    for (i in 0..n) {
        diff.data[i] *= scale
    }
    gemm(NoTrans, NoTrans, 1.0, a, x, -1.0, diff)

    let xnorm = dlange(MaxColumnSum, n, 1, x.data, 1)
    let resid = dlange(MaxColumnSum, n, 1, diff.data, 1) / xnorm
    if (resid > tol || resid.isNaN()) {
        throw Exception("${name}: unexpected result; resid=${resid}, want<=${tol}")
    }
    return true
}

func copy(dst: Array<Int>, src: Array<Int>): Int64 {
    let minSize = min(src.size, dst.size)
    dst[0..minSize] = src[0..minSize]
    return minSize
}
