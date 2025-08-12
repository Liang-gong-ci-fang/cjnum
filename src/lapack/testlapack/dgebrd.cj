package cjnum.lapack.testlapack

import cjnum.lapack.*
import cjnum.floats.*
import std.math.*
import std.random.*

class dgebrdTestStruct {
    dgebrdTestStruct(
        let m: Int64,
        let n: Int64,
        let lda: Int64
    ) {}
}

public func dgebrdTest(lFloat64: LFloat64): Bool {
    let rnd = Random(1)
    for (test in [
            dgebrdTestStruct(100, 100, 0),
            dgebrdTestStruct(100, 150, 0),
            dgebrdTestStruct(150, 100, 0),
            dgebrdTestStruct(100, 100, 200),
            dgebrdTestStruct(100, 150, 200),
            dgebrdTestStruct(150, 100, 200),
            dgebrdTestStruct(300, 300, 0),
            dgebrdTestStruct(300, 400, 0),
            dgebrdTestStruct(400, 300, 0),
            dgebrdTestStruct(300, 300, 500),
            dgebrdTestStruct(300, 400, 500),
            dgebrdTestStruct(300, 400, 500)
        ]) {
        let m = test.m
        let n = test.n
        var lda = test.lda
        if (lda == 0) {
            lda = n
        }
        let minmn = min(m, n)
        let a = Array<Float64>(m * lda, repeat: 0.0)

        for (i in 0..a.size - 1) {
            a[i] = rnd.nextGaussianFloat64()
        }

        let d = Array<Float64>(minmn, repeat: 0.0)
        let e = Array<Float64>(minmn - 1, repeat: 0.0)
        let tauP = Array<Float64>(minmn, repeat: 0.0)
        let tauQ = Array<Float64>(minmn, repeat: 0.0)
        var work = Array<Float64>(max(m, n), repeat: 0.0)
        for (i in 0..work.size - 1) {
            work[i] = Float64.NaN
        }

        let aCopy = Array<Float64>(a.size, repeat: 0.0)
        copy(aCopy, a)

        lFloat64.dgebd2(m, n, a, lda, d, e, tauQ, tauP, work)
        let aAns = Array<Float64>(a.size, repeat: 0.0)
        copy(aAns, a)
        let dAns = Array<Float64>(d.size, repeat: 0.0)
        copy(dAns, d)
        let eAns = Array<Float64>(e.size, repeat: 0.0)
        copy(eAns, e)
        let tauQAns = Array<Float64>(tauQ.size, repeat: 0.0)
        copy(tauQAns, tauQ)
        let tauPAns = Array<Float64>(tauP.size, repeat: 0.0)
        copy(tauPAns, tauP)

        var lwork = -1
        copy(a, aCopy)
        lFloat64.dgebrd(m, n, a, lda, d, e, tauQ, tauP, work, lwork)
        work = Array<Float64>(Int64(work[0]), repeat: 0.0)
        lwork = work.size
        for (i in 0..work.size) {
            work[i] = Float64.NaN
        }
        for (i in 0..d.size) {
            d[i] = Float64.NaN
        }
        for (i in 0..e.size) {
            e[i] = Float64.NaN
        }
        for (i in 0..tauQ.size) {
            tauQ[i] = Float64.NaN
        }
        for (i in 0..tauP.size) {
            tauP[i] = Float64.NaN
        }
        lFloat64.dgebrd(m, n, a, lda, d, e, tauQ, tauP, work, lwork)

        if (!fEqualApprox(a, aAns, 1e-10)) {
            throw Exception("a mismatch")
        }
        if (!fEqualApprox(d, dAns, 1e-10)) { 
            throw Exception("d mismatch")
        }
        if (!fEqualApprox(e, eAns, 1e-10)) {
            throw Exception("e mismatch")
        }
        if (!fEqualApprox(tauQ, tauQAns, 1e-10)) {
            throw Exception("tauQ mismatch")
        }
        if (!fEqualApprox(tauP, tauPAns, 1e-10)) {
            throw Exception("tauP mismatch")
        }

        lwork--
        copy(a, aCopy)
        for (i in 0..d.size) {
            d[i] = 0.0
        }
        for (i in 0..e.size) {
            e[i] = 0.0
        }
        for (i in 0..tauP.size) {
            tauP[i] = 0.0
        }
        for (i in 0..tauQ.size) {
            tauQ[i] = 0.0
        }
        lFloat64.dgebrd(m, n, a, lda, d, e, tauQ, tauP, work, lwork)

        if (!fEqualApprox(a, aAns, 1e-10)) {
            throw Exception("a mismatch")
        }
        if (!fEqualApprox(d, dAns, 1e-10)) {
            throw Exception("d mismatch")
        }
        if (!fEqualApprox(e, eAns, 1e-10)) {
            throw Exception("e mismatch")
        }
        if (!fEqualApprox(tauQ, tauQAns, 1e-10)) {
            throw Exception("tauQ mismatch")
        }
        if (!fEqualApprox(tauP, tauPAns, 1e-10)) {
            throw Exception("tauP mismatch")
        }
    }
    return true
}
