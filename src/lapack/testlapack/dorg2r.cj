package cjnum.lapack.testlapack

import cjnum.floats.scalar.*
import cjnum.lapack.*
import std.math.*
import std.random.*

class dorg2rTestStruct {
    dorg2rTestStruct(
        let m: Int64,
        let n: Int64,
        let k: Int64,
        let lda: Int64
    ) {}
}

public func dorg2rTest(lFloat64: LFloat64): Bool {
    let rnd = Random(1)
    var ti = 0
    for (test in [
            dorg2rTestStruct(3, 3, 0, 0),
            dorg2rTestStruct(4, 3, 0, 0),
            dorg2rTestStruct(3, 3, 2, 0),
            dorg2rTestStruct(4, 3, 2, 0),
            dorg2rTestStruct(5, 5, 0, 20),
            dorg2rTestStruct(5, 5, 3, 20),
            dorg2rTestStruct(10, 5, 0, 20),
            dorg2rTestStruct(10, 5, 2, 20)
        ]) {
        var m = test.m
        var n = test.n
        var lda = test.lda
        if (lda == 0) {
            lda = test.n
        }

        let a = Array<Float64>(m * lda, repeat: 0.0)
        for (i in 0..a.size) {
            a[i] = rnd.nextGaussianFloat64()
        }

        let tau = Array<Float64>(min(m, n), repeat: 0.0)
        var work = Array<Float64>(1, repeat: 0.0)
        lFloat64.dgeqrf(m, n, a, lda, tau, work, -1)
        work = Array<Float64>(Int64(work[0]), repeat: 0.0)
        lFloat64.dgeqrf(m, n, a, lda, tau, work, work.size)

        var k = test.k
        if (k == 0) {
            k = n
        }
        let q = constructQK("QR", m, n, k, a, lda, tau)

        lFloat64.dorg2r(m, n, k, a, lda, tau[..k], work)

        var same = true
        var loop = true

        for (i in 0..m) {
            if (!loop) {
                break
            }
            for (j in 0..n) {
                if (!equalWithinAbsOrRel(q.data[i * q.stride + j], a[i * lda + j], 1e-12, 1e-12)) {
                    same = false
                    loop = false
                    break
                }
            }
        }
        if (!same) {
            throw Exception("Case ${ti}: Q mismatch")
        }
        ti++
    }
    return true
}
