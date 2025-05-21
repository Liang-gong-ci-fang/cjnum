package cjnum.lapack.testlapack

import std.math.*
import std.random.*
import cjnum.lapack.*

class dgebd2TestStruct {
    let m: Int
    let n: Int
    let lda: Int

    init(m: Int, n: Int, lda: Int) {
        this.m = m
        this.n = n
        this.lda = lda
    }
}

public func dgebd2Test(lFloat64: LFloat64): Bool {
    let rnd = Random(1)
    for (test in [
            dgebd2TestStruct(3, 4, 0),
            dgebd2TestStruct(4, 3, 0),
            dgebd2TestStruct(3, 4, 10),
            dgebd2TestStruct(4, 3, 10)
        ]) {
        let m = test.m
        let n = test.n
        var lda = test.lda
        if (lda == 0) {
            lda = n
        }

        let a = Array<Float64>(m * lda, item: 0.0)
        for (i in 0..a.size) {
            a[i] = rnd.nextGaussianFloat64()
        }

        let aCopy = Array<Float64>(a.size, item: 0.0)
        copy(aCopy, a)

        let nb = min(m, n)
        let d = nanSlice(nb)
        let e = nanSlice(nb - 1)

        let tauP = nanSlice(nb)
        let tauQ = nanSlice(nb)

        let work = nanSlice(max(m, n))

        lFloat64.dgebd2(m, n, a, lda, d, e, tauQ, tauP, work)

        checkBidiagonal(lFloat64, m, n, nb, a, lda, d, e, tauP, tauQ, aCopy)
    }
    return true
}
