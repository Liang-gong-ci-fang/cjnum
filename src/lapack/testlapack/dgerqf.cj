package cjnum.lapack.testlapack

import cjnum.blas.blas64.*
import cjnum.blas.*
import cjnum.lapack.*
import std.math.*
import std.random.*

class dgerqfTestStruct {
    dgerqfTestStruct(
        let name!: String,
        let length!: Int64
    ) {}
}

public func dgerqfTest(lFloat64: LFloat64): Bool {
    let rnd = Random(1)
    for (m in [0, 1, 2, 3, 4, 5, 6, 12, 129, 160]) {
        for (n in [0, 1, 2, 3, 4, 5, 6, 12, 129, 160]) {
            for (lda in [max(1, n), n + 4]) {
                if (!dgerqfTest(lFloat64, rnd, m, n, lda)) {
                    return false
                }
            }
        }
    }
    return true
}

func dgerqfTest(lFloat64: LFloat64, rnd: Random, m: Int64, n: Int64, lda: Int64): Bool {
    const tol = 1e-14

    let a = randomGeneral(m, n, lda, rnd)
    let aCopy = cloneGeneral(a)

    let k = min(m, n)
    let tau = Array<Float64>(k, {_ => rnd.nextFloat64()})

    var work = [0.0]
    lFloat64.dgerqf(m, n, a.data, a.stride, tau, work, -1)
    let lwkopt = Int64(work[0])
    for (wk in [
            dgerqfTestStruct(name: "short", length: m),
            dgerqfTestStruct(name: "medium", length: lwkopt - 1),
            dgerqfTestStruct(name: "long", length: lwkopt)
        ]) {
        let name = "m=${m},n=${n},lda=${lda},work=${wk.name}"

        let lwork = max(max(1, m), wk.length)
        work = Array<Float64>(lwork, {_ => rnd.nextFloat64()})

        copyGeneral(a, aCopy)
        lFloat64.dgerqf(m, n, a.data, a.stride, tau, work, lwork)

        let q = constructQ("RQ", m, n, a.data, a.stride, tau)

        var resid = residualOrthogonal(q, false)
        if (resid > tol * Float64(m)) {
            throw Exception("Case ${name}: Q not orthogonal; resid=${resid}, want<=${tol*Float64(m)}")
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
    }
    return true
}
