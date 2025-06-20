package cjnum.lapack.testlapack

import cjnum.blas.blas64.*
import cjnum.blas.*
import cjnum.floats.*
import cjnum.lapack.*
import std.math.*
import std.random.*

public func dgerq2Test(lFloat64: LFloat64): Bool {
    let rnd = Random(1)
    for (m in [0, 1, 2, 3, 4, 5, 6, 12, 23]) {
        for (n in [0, 1, 2, 3, 4, 5, 6, 12, 23]) {
            for (lda in [max(1, n), n + 4]) {
                if (!dgerq2Test(lFloat64, rnd, m, n, lda)) {
                    return false
                }
            }
        }
    }
    return true
}

func dgerq2Test(lFloat64: LFloat64, rnd: Random, m: Int64, n: Int64, lda: Int64): Bool {
    const tol = 1e-14

    let name = "m=${m},n=${n},lda=${lda}"

    let a = randomGeneral(m, n, lda, rnd)
    let aCopy = cloneGeneral(a)

    let k = min(m, n)
    let tau = Array<Float64>(k, {_ => rnd.nextFloat64()})

    let work = Array<Float64>(m, {_ => rnd.nextFloat64()})

    lFloat64.dgerq2(m, n, a.data, a.stride, tau, work)

    let q = constructQ("RQ", m, n, a.data, a.stride, tau)

    var resid = residualOrthogonal(q, false)
    if (resid > tol) {
        throw Exception("Case ${name}: Q not orthogonal; resid=${resid}, want<=${tol}")
    }

    let r = zeros(m, n, n)
    for (i in 0..m) {
        let off = m - n
        for (j in max(0, i - off)..n) {
            r.data[i * r.stride + j] = a.data[i * a.stride + j]
        }
    }
    let qra = cloneGeneral(aCopy)
    gemm(NoTrans, NoTrans, 1.0, r, q, -1.0, qra)
    resid = dlange(MaxColumnSum, qra.rows, qra.cols, qra.data, qra.stride)
    if (resid > tol * Float64(m)) {
        throw Exception("Case ${name}: |R*Q - A|=${resid}, want<=${tol*Float64(m)}")
    }
    return true
}
