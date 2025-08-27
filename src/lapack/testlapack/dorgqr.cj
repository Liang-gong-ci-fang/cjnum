package cjnum.lapack.testlapack

import cjnum.floats.*
import cjnum.lapack.*
import std.random.*

class dorgqrTestStruct {
    dorgqrTestStruct(
        let m: Int64,
        let n: Int64,
        let k: Int64,
        let lda: Int64
    ) {}
}

public func dorgqrTest(lFloat64: LFloat64): Bool {
    let rnd = Random(1)

    for (test in [
            dorgqrTestStruct(10, 10, 10, 0),
            dorgqrTestStruct(10, 10, 10, 20),
            dorgqrTestStruct(30, 10, 10, 0),
            dorgqrTestStruct(30, 20, 10, 20),
            dorgqrTestStruct(100, 100, 100, 0),
            dorgqrTestStruct(100, 100, 50, 0),
            dorgqrTestStruct(130, 100, 100, 0),
            dorgqrTestStruct(130, 100, 50, 0),
            dorgqrTestStruct(100, 100, 100, 150),
            dorgqrTestStruct(100, 100, 50, 150),
            dorgqrTestStruct(130, 100, 100, 150),
            dorgqrTestStruct(130, 100, 50, 150),
            dorgqrTestStruct(200, 200, 200, 0),
            dorgqrTestStruct(200, 200, 150, 0),
            dorgqrTestStruct(230, 200, 200, 0),
            dorgqrTestStruct(230, 200, 150, 0),
            dorgqrTestStruct(200, 200, 200, 250),
            dorgqrTestStruct(200, 200, 150, 250),
            dorgqrTestStruct(230, 200, 200, 250),
            dorgqrTestStruct(230, 200, 150, 250)
        ]) {
        let m = test.m
        let n = test.n
        let k = test.k
        var lda = test.lda
        if (lda == 0) {
            lda = n
        }
        let a = Array<Float64>(
            m * lda,
            {
                _ => rnd.nextFloat64()
            }
        )
        var work = Array<Float64>(1, repeat: 0.0)
        let tau = Array<Float64>(
            min(m, n),
            {
                _ => Float64.NaN
            }
        )

        lFloat64.dgeqrf(m, n, a, lda, tau, work, -1)
        work = Array<Float64>(Int64(work[0]), repeat: 0.0)
        lFloat64.dgeqrf(m, n, a, lda, tau, work, work.size)

        let aUnblocked = Array<Float64>(a.size, repeat: 0.0)
        copy(aUnblocked, a)
        for (i in 0..work.size) {
            work[i] = Float64.NaN
        }
        lFloat64.dorg2r(m, n, k, aUnblocked, lda, tau[..k], work)

        for (i in 0..work.size) {
            work[i] = Float64.NaN
        }
        lFloat64.dorgqr(m, n, k, a, lda, tau[..k], work, work.size)
        if (!fEqualApprox(a, aUnblocked, 1e-10)) {
            throw Exception("Q Mismatch. m = ${m}, n = ${n}, k = ${k}, lda = ${lda}")
        }
    }
    return true
}
