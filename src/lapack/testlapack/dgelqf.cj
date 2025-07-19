package cjnum.lapack.testlapack

import cjnum.floats.*
import cjnum.lapack.*
import std.random.*

class dgelqfTestStruct {
    dgelqfTestStruct(
        var m: Int64,
        var n: Int64,
        var lda: Int64
    ) {}
}

public func dgelqfTest(lFloat64: LFloat64): Bool {
    const tol = 1e-12
    let rnd = Random(1)
    var c = 0
    for (test in [
            dgelqfTestStruct(10, 5, 0),
            dgelqfTestStruct(5, 10, 0),
            dgelqfTestStruct(10, 10, 0),
            dgelqfTestStruct(300, 5, 0),
            dgelqfTestStruct(3, 500, 0),
            dgelqfTestStruct(200, 200, 0),
            dgelqfTestStruct(300, 200, 0),
            dgelqfTestStruct(204, 300, 0),
            dgelqfTestStruct(1, 3000, 0),
            dgelqfTestStruct(3000, 1, 0),
            dgelqfTestStruct(10, 5, 30),
            dgelqfTestStruct(5, 10, 30),
            dgelqfTestStruct(10, 10, 30),
            dgelqfTestStruct(300, 5, 500),
            dgelqfTestStruct(3, 500, 600),
            dgelqfTestStruct(200, 200, 300),
            dgelqfTestStruct(300, 200, 300),
            dgelqfTestStruct(204, 300, 400),
            dgelqfTestStruct(1, 3000, 4000),
            dgelqfTestStruct(3000, 1, 4000)
        ]) {
        let m = test.m
        let n = test.n
        var lda = test.lda
        if (lda == 0) {
            lda = n
        }
        let a = Array<Float64>(m * lda, {_ => rnd.nextGaussianFloat64()})

        let aCopy = Array<Float64>(a.size, item: 0.0)
        copy(aCopy, a)

        let tau = Array<Float64>(n, {_ => rnd.nextGaussianFloat64()})

        let want = Array<Float64>(a.size, item: 0.0)
        copy(want, a)
        lFloat64.dgelq2(m, n, want, lda, tau, Array<Float64>(m, item: 0.0))

        for (wl in [minimumWork, mediumWork, optimumWork]) {
            copy(a, aCopy)

            var lwork = 0
            if (wl == minimumWork) {
                lwork = m
            } else if (wl == mediumWork) {
                let work = Array<Float64>(1, item: 0.0)
                lFloat64.dgelqf(m, n, a, lda, tau, work, -1)
                lwork = Int64(work[0]) - 2 * m
            } else if (wl == optimumWork) {
                let work = Array<Float64>(1, item: 0.0)
                lFloat64.dgelqf(m, n, a, lda, tau, work, -1)
                lwork = Int64(work[0])
            }
            let work = Array<Float64>(lwork, item: 0.0)

            lFloat64.dgelqf(m, n, a, lda, tau, work, work.size)
            if (!fEqualApprox(want, a, tol)) {
                throw Exception("Case ${c}, workspace type ${wl}, unexpected result")
            }
        }
        c++
    }
    return true
}
