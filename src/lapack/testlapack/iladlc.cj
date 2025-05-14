package cjnum.lapack.testlapack

import cjnum.lapack.*

class IladlcTestStruct {
    IladlcTestStruct(
        let a!: Array<Float64>,
        let m!: Int64,
        let n!: Int64,
        let lda!: Int64,
        let ans!: Int64
    ) {}
}

public func iladlcTest(lFloat64: LFloat64): Bool {
    let tests = [
        IladlcTestStruct(
            a: [0.0, 0.0, 0.0, 0.0],
            m: 1,
            n: 1,
            lda: 2,
            ans: -1,
        ),
        IladlcTestStruct(
            a: [0.0, 0.0, 0.0, 0.0],
            m: 2,
            n: 2,
            lda: 2,
            ans: -1
        ),
        IladlcTestStruct(
            a: [0.0, 0.0, 0.0, 0.0],
            m: 4,
            n: 1,
            lda: 1,
            ans: -1
        ),
        IladlcTestStruct(
            a: [0.0, 0.0, 0.0, 0.0],
            m: 1,
            n: 4,
            lda: 4,
            ans: -1
        ),
        IladlcTestStruct(
            a: [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0],
            m: 2,
            n: 4,
            lda: 4,
            ans: 3
        ),
        IladlcTestStruct(
            a: [1.0, 2.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            m: 2,
            n: 4,
            lda: 4,
            ans: 2
        ),
        IladlcTestStruct(
            a: [0.0, 0.0, 3.0, 4.0, 0.0, 0.0, 0.0, 0.0],
            m: 2,
            n: 2,
            lda: 4,
            ans: -1
        )
    ]
    for (i in 0..tests.size) {
        let test = tests[i]
        let ans = lFloat64.iladlc(test.m, test.n, test.a, test.lda)
        if (ans != test.ans) {
            throw Exception("Column mismatch case ${i}. Want: ${test.ans}, got: ${ans}")
        }
    }
    return true
}