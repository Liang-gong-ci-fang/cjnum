package cjnum.lapack.testlapack

import std.math.*
import std.random.*
import cjnum.lapack.*

class dlabrdTestStruct {
    dlabrdTestStruct(
        let m: Int,
        let n: Int,
        let nb: Int,
        let lda: Int,
        let ldx: Int,
        let ldy: Int
    ) {}
}

public func dlabrdTest(lFloat64: LFloat64): Bool {
    let rnd = Random(1)
    for (test in [
            dlabrdTestStruct(4, 5, 2, 0, 0, 0),
            dlabrdTestStruct(4, 5, 4, 0, 0, 0),
            dlabrdTestStruct(5, 5, 2, 0, 0, 0),
            dlabrdTestStruct(5, 5, 5, 0, 0, 0),
            dlabrdTestStruct(5, 4, 2, 0, 0, 0),
            dlabrdTestStruct(5, 4, 4, 0, 0, 0),
            dlabrdTestStruct(4, 5, 2, 10, 11, 12),
            dlabrdTestStruct(4, 5, 4, 10, 11, 12),
            dlabrdTestStruct(5, 5, 2, 10, 11, 12),
            dlabrdTestStruct(5, 5, 5, 10, 11, 12),
            dlabrdTestStruct(5, 4, 2, 10, 11, 12),
            dlabrdTestStruct(5, 4, 4, 10, 11, 12),
            dlabrdTestStruct(4, 5, 2, 11, 12, 10),
            dlabrdTestStruct(4, 5, 4, 11, 12, 10),
            dlabrdTestStruct(5, 5, 2, 11, 12, 10),
            dlabrdTestStruct(5, 5, 5, 11, 12, 10),
            dlabrdTestStruct(5, 4, 2, 11, 12, 10),
            dlabrdTestStruct(5, 4, 4, 11, 12, 10),
            dlabrdTestStruct(4, 5, 2, 12, 11, 10),
            dlabrdTestStruct(4, 5, 4, 12, 11, 10),
            dlabrdTestStruct(5, 5, 2, 12, 11, 10),
            dlabrdTestStruct(5, 5, 5, 12, 11, 10),
            dlabrdTestStruct(5, 4, 2, 12, 11, 10),
            dlabrdTestStruct(5, 4, 4, 12, 11, 10)
        ]) {
        let m = test.m
        let n = test.n
        let nb = test.nb
        var lda = test.lda
        if (lda == 0) {
            lda = n
        }
        var ldy = test.ldy
        if (ldy == 0) {
            ldy = nb
        }
        var ldx = test.ldx
        if (ldx == 0) {
            ldx = nb
        }
        let a = Array<Float64>(m * lda, {_ => rnd.nextGaussianFloat64()})

        let d = Array<Float64>(nb, {_ => Float64.NaN})

        let e = Array<Float64>(nb, {_ => Float64.NaN})

        let tauP = Array<Float64>(nb, {_ => Float64.NaN})

        let tauQ = Array<Float64>(nb, {_ => Float64.NaN})

        let x = Array<Float64>(m * ldx, {_ => rnd.nextGaussianFloat64()})

        let y = Array<Float64>(n * ldy, {_ => rnd.nextGaussianFloat64()})

        let aCopy = Array<Float64>(a.size, item: 0.0)
        copy(aCopy, a)

        lFloat64.dlabrd(m, n, nb, a, lda, d, e, tauQ, tauP, x, ldx, y, ldy)

        if (m >= n && nb == n) {
            tauP[n - 1] = 0.0
        }
        if (m < n && nb == m) {
            tauQ[m - 1] = 0.0
        }
        checkBidiagonal(m, n, nb, a, lda, d, e, tauP, tauQ, aCopy) //TODO 在此处中精度存在问题，误差为10
    }
    return true
}
