package cjnum.lapack.testlapack

import cjnum.floats.*
import cjnum.lapack.*
import std.random.*

class dgeqrfTestStruct {
    dgeqrfTestStruct(
        let m: Int64,
        let n: Int64,
        let lda: Int64
    ) {}
}

public func dgeqrfTest(lFloat64: LFloat64): Bool {
    const tol = 1e-12
    let rnd = Random(1)
    var c = 0
    for (test in [
            dgeqrfTestStruct(10, 5, 0),
            dgeqrfTestStruct(5, 10, 0),
            dgeqrfTestStruct(10, 10, 0),
            dgeqrfTestStruct(300, 5, 0),
            dgeqrfTestStruct(3, 500, 0),
            dgeqrfTestStruct(200, 200, 0),
            dgeqrfTestStruct(300, 200, 0),
            dgeqrfTestStruct(204, 300, 0),
            dgeqrfTestStruct(1, 3000, 0),
            dgeqrfTestStruct(3000, 1, 0),
            dgeqrfTestStruct(10, 5, 20),
            dgeqrfTestStruct(5, 10, 20),
            dgeqrfTestStruct(10, 10, 20),
            dgeqrfTestStruct(300, 5, 400),
            dgeqrfTestStruct(3, 500, 600),
            dgeqrfTestStruct(200, 200, 300),
            dgeqrfTestStruct(300, 200, 300),
            dgeqrfTestStruct(204, 300, 400),
            dgeqrfTestStruct(1, 3000, 4000),
            dgeqrfTestStruct(3000, 1, 4000)
        ]) {
        let m = test.m
        let n = test.n
        var lda = test.lda
        if (lda == 0) {
            lda = test.n
        }

        let a = Array<Float64>(m * lda,  { _ => rnd.nextGaussianFloat64()})

        let aCopy = Array<Float64>(a.size, repeat: 0.0)
        copy(aCopy, a)

        let k = min(m, n)
        let tau = Array<Float64>(k, { _ => rnd.nextFloat64()})

        let want = Array<Float64>(a.size, repeat: 0.0)
        copy(want, a)
        lFloat64.dgeqr2(m, n, want, lda, tau, Array<Float64>(n, repeat: 0.0))

        for (wl in [minimumWork, mediumWork, optimumWork]) {
            copy(a, aCopy)

            var lwork: Int64 = 0
            if (wl == minimumWork) {
                lwork = n
            } else if (wl == mediumWork) {
                let work = [0.0]
                lFloat64.dgeqrf(m, n, a, lda, tau, work, -1)
                lwork = Int64(work[0]) - 2 * n
            } else if (wl == optimumWork) {
                let work = [0.0]
                lFloat64.dgeqrf(m, n, a, lda, tau, work, -1)
                lwork = Int64(work[0])
            }
            let work = Array<Float64>(lwork, repeat: 0.0)

            lFloat64.dgeqrf(m, n, a, lda, tau, work, work.size)

            if (!fEqualApprox(want, a, tol)) {
                throw Exception("Case ${c}, workspace ${wl}, unexpected result.")
            }
        }
        c++
    }
    return true
}
