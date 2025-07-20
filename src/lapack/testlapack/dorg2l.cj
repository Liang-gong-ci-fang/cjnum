package cjnum.lapack.testlapack

import cjnum.blas.blas64.General as General64
import cjnum.blas.blas64.*
import cjnum.lapack.*
import std.math.*
import std.random.*

class dorg2lTestStruct {
    dorg2lTestStruct(
        let m: Int64,
        let n: Int64,
        let k: Int64,
        let lda: Int64
    ) {}
}

public func dorg2lTest(lFloat64: LFloat64): Bool {
    const tol = 1e-14

    let rnd = Random(1)
    for (test in [
            dorg2lTestStruct(5, 4, 3, 0),
            dorg2lTestStruct(5, 4, 4, 0),
            dorg2lTestStruct(3, 3, 2, 0),
            dorg2lTestStruct(5, 5, 5, 0),
            dorg2lTestStruct(5, 4, 3, 11),
            dorg2lTestStruct(5, 4, 4, 11),
            dorg2lTestStruct(3, 3, 2, 11),
            dorg2lTestStruct(5, 5, 5, 11)
        ]) {
        let m = test.m
        let n = test.n
        let k = test.k
        var lda = test.lda
        if (lda == 0) {
            lda = n
        }

        let a = Array<Float64>(m * lda, repeat: 0.0)
        for (i in 0..a.size) {
            a[i] = rnd.nextGaussianFloat64()
        }
        let tau = nanSlice(max(m, n))
        let work = Array<Float64>(n, repeat: 0.0)
        lFloat64.dgeql2(m, n, a, lda, tau, work)

        lFloat64.dorg2l(m, n, k, a, lda, tau[n - k..], work)

        let q = General64(rows: m, cols: n, data: a, stride: lda)
        if (residualOrthogonal(q, false) > tol) {
            throw Exception(
                "Case m=${m}, n=${n}, k=${k}, lda=${lda}: columns of Q not orthonormal; resid=${residualOrthogonal(q, false)}, want<=${tol}"
            )
        }
    }
    return true
}
