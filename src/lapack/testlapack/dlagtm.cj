package cjnum.lapack.testlapack

import cjnum.blas.blas64.*
import cjnum.blas.*
import cjnum.floats.*
import cjnum.lapack.*
import std.math.*
import std.random.*

public func dlagtmTest(lFloat64: LFloat64): Bool {
    let rnd = Random(1)
    for (trans in [NoTrans, Trans, ConjTrans]) {
        for (m in [0, 1, 2, 3, 4, 5, 10]) {
            for (n in [0, 1, 2, 3, 4, 5, 10]) {
                for (ldb in [max(1, n), n + 3]) {
                    for (ldc in [max(1, n), n + 4]) {
                        for (alpha in [0.0, 1.0, rnd.nextGaussianFloat64()]) {
                            for (beta in [0.0, 1.0, rnd.nextGaussianFloat64()]) {
                                if (!dlagtmTest(lFloat64, rnd, trans, m, n, ldb, ldc, alpha, beta)) {
                                    return false
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return true
}

func dlagtmTest(lFloat64: LFloat64, rnd: Random, trans: Transpose, m: Int64, n: Int64, ldb: Int64, ldc: Int64,
    alpha: Float64, beta: Float64): Bool {
    const tol = 1e-14
    const extra = 10

    let name = "Case m=${m},n=${n},ldb=${ldb},ldc=${ldc},alpha=${alpha},beta=${beta}"

    let dl = randomSlice(n + extra, rnd)
    let dlCopy = Array<Float64>(dl.size, repeat: 0.0)
    copy(dlCopy, dl)

    let d = randomSlice(n + 1 + extra, rnd)
    let dCopy = Array<Float64>(d.size, repeat: 0.0)
    copy(dCopy, d)

    let du = randomSlice(n + extra, rnd)
    let duCopy = Array<Float64>(du.size, repeat: 0.0)
    copy(duCopy, du)

    let b = randomGeneral(m, n, ldb, rnd)
    let bCopy = cloneGeneral(b)

    let got = randomGeneral(m, n, ldc, rnd)
    let want = cloneGeneral(got)

    if (m == 0 || n == 0) {
        lFloat64.dlagtm(trans, m, n, alpha, dl, d, du, b.data, b.stride, beta, got.data, got.stride)
        if (!fSame(dl, dlCopy)) {
            throw Exception("${name}: unexpected modification in dl")
        }
        if (!fSame(d, dCopy)) {
            throw Exception("${name}: unexpected modification in d")
        }
        if (!fSame(du, duCopy)) {
            throw Exception("${name}: unexpected modification in du")
        }
        if (!fSame(b.data, bCopy.data)) {
            throw Exception("${name}: unexpected modification in B")
        }
        if (!fSame(got.data, want.data)) {
            throw Exception("${name}: unexpected modification in C")
        }
        return true
    }

    lFloat64.dlagtm(trans, m, n, alpha, dl, d, du, b.data, b.stride, beta, got.data, got.stride)

    if (!fSame(dl, dlCopy)) {
        throw Exception("${name}: unexpected modification in dl")
    }
    if (!fSame(d, dCopy)) {
        throw Exception("${name}: unexpected modification in d")
    }
    if (!fSame(du, duCopy)) {
        throw Exception("${name}: unexpected modification in du")
    }
    if (!fSame(b.data, bCopy.data)) {
        throw Exception("${name}: unexpected modification in B")
    }

    let a = zeros(m, m, m)
    for (i in 0..m - 1) {
        a.data[i * a.stride + i] = d[i]
        a.data[i * a.stride + i + 1] = du[i]
        a.data[(i + 1) * a.stride + i] = dl[i]
    }
    a.data[(m - 1) * a.stride + m - 1] = d[m - 1]

    gemm(trans, NoTrans, alpha, a, b, beta, want)

    for (i in 0..m) {
        for (j in 0..n) {
            got.data[i * got.stride + j] -= want.data[i * want.stride + j]
        }
    }
    let diff = dlange(MaxColumnSum, got.rows, got.cols, got.data, got.stride)
    if (diff > tol) {
        throw Exception("${name}: unexpected result; diff=${diff}, want<=${tol}")
    }
    return true
}
