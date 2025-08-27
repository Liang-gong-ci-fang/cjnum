package cjnum.lapack.testlapack

import cjnum.blas.blas64.*
import cjnum.blas.*
import cjnum.lapack.*
import std.random.*

public func dgetc2Test(lFloat64: LFloat64): Bool {
    let rnd = Random(1)
    for (n in [0, 1, 2, 3, 4, 5, 10, 20]) {
        for (lda in [n, n + 5]) {
            if (!dgetc2Test(lFloat64, rnd, n, lda, false)) {
                return false
            }
            if (!dgetc2Test(lFloat64, rnd, n, lda, true)) {
                return false
            }
        }
    }
    return true
}

func dgetc2Test(lFloat64: LFloat64, rnd: Random, n: Int64, lda: Int64, perturb: Bool): Bool {
    const tol = 1e-14

    let name = "n=${n},lda=${lda},perturb=${perturb}"

    var l = randomGeneral(n, n, max(1, n), rnd)
    for (i in 0..n) {
        l.data[i * l.stride + i] = 1.0
        for (j in i + 1..n) {
            l.data[i * l.stride + j] = 0.0
        }
    }

    var u = randomGeneral(n, n, max(1, n), rnd)
    for (i in 0..n) {
        for (j in 0..i) {
            u.data[i * u.stride + j] = 0.0
        }
    }
    if (perturb && n > 0) {
        let i = rnd.nextInt64(n)
        u.data[i * u.stride + i] = 0.0
    }

    let a = zeros(n, n, max(1, lda))
    gemm(NoTrans, NoTrans, 1.0, l, u, 0.0, a)

    let ipiv = Array<Int64>(n, repeat: 0)
    let jpiv = Array<Int64>(n, repeat: 0)
    for (i in 0..n) {
        ipiv[i] = -1
        jpiv[i] = -1
    }

    let lu = cloneGeneral(a)
    let k = lFloat64.dgetc2(n, lu.data, lu.stride, ipiv, jpiv)

    if (n == 0) {
        return true
    }

    if (perturb && k < 0) {
        throw Exception("${name}: expected matrix perturbation")
    }

    for (i in 0..n) {
        if (ipiv[i] < 0) {
            throw Exception("${name}: ipiv[${i}] is not set")
        }
        if (jpiv[i] < 0) {
            throw Exception("${name}: jpiv[${i}] is not set")
        }
    }

    l = zeros(n, n, n)
    u = zeros(n, n, n)
    for (i in 0..n) {
        for (j in 0..i) {
            l.data[i * l.stride + j] = lu.data[i * lu.stride + j]
        }
        l.data[i * l.stride + i] = 1.0
        for (j in i..n) {
            u.data[i * u.stride + j] = lu.data[i * lu.stride + j]
        }
    }
    let diff = zeros(n, n, n)
    gemm(NoTrans, NoTrans, 1.0, l, u, 0.0, diff)

    for (i in n - 1..=0: -1) {
        let ipv = ipiv[i]
        if (ipv != i) {
            let row1 = Vector(n: n, data: diff.data[i * diff.stride..], inc: 1)
            let row2 = Vector(n: n, data: diff.data[ipv * diff.stride..], inc: 1)
            swap(row1, row2)
        }
        let jpv = jpiv[i]
        if (jpv != i) {
            let col1 = Vector(n: n, data: diff.data[i..], inc: diff.stride)
            let col2 = Vector(n: n, data: diff.data[jpv..], inc: diff.stride)
            swap(col1, col2)
        }
    }

    for (i in 0..n) {
        for (j in 0..n) {
            diff.data[i * diff.stride + j] -= a.data[i * a.stride + j]
        }
    }
    let resid = dlange(MaxColumnSum, n, n, diff.data, diff.stride)
    if (resid > tol || resid.isNaN()) {
        throw Exception("${name}: unexpected result; resid=${resid}, want<=${tol}")
    }
    return true
}
