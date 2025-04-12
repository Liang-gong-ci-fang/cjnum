package cjnum.blas.testblas

import cjnum.blas.*
import cjnum.floats.*

class dsyrTestStruct {
    dsyrTestStruct(
        let ul!: Uplo,
        let n!: Int64,
        let a!: Array<Array<Float64>>,
        let x!: Array<Float64>,
        let alpha!: Float64,
        let ans!: Array<Array<Float64>>
    ) {}
}

public func dsyrTest(nFloat64: NFloat64): Bool {
    let tests = [
        dsyrTestStruct(
            ul: Upper,
            n: 4,
            a: [
                [10.0, 2.0, 0.0, 1.0],
                [0.0, 1.0, 2.0, 3.0],
                [0.0, 0.0, 9.0, 15.0],
                [0.0, 0.0, 0.0, -6.0]
            ],
            x: [1.0, 2.0, 0.0, 5.0],
            alpha: 8.0,
            ans: [
                [18.0, 18.0, 0.0, 41.0],
                [0.0, 33.0, 2.0, 83.0],
                [0.0, 0.0, 9.0, 15.0],
                [0.0, 0.0, 0.0, 194.0]
            ]
        ),
        dsyrTestStruct(
            ul: Lower,
            n: 3,
            a: [
                [10.0, 2.0, 0.0],
                [4.0, 1.0, 2.0],
                [2.0, 7.0, 9.0]
            ],
            x: [3.0, 0.0, 5.0],
            alpha: 8.0,
            ans: [
                [82.0, 2.0, 0.0],
                [4.0, 1.0, 2.0],
                [122.0, 7.0, 209.0]
            ]
        )
    ]

    for (i in 0..tests.size) {
        let incTest = { incX: Int64, extra: Int64 =>
            let xnew = makeIncremented(tests[i].x, incX, extra)
            let aFlat = flatten(tests[i].a)
            let ans = flatten(tests[i].ans)
            let lda = tests[i].n
            nFloat64.dsyr(tests[i].ul, tests[i].n, tests[i].alpha, xnew, incX, aFlat, lda)
            if (!dSliceTolEqual(aFlat, ans)) {
                throw Exception("Case ${i}, idx ${incX}: Want ${ans}, got ${aFlat}.")
            }
        }

		incTest(1, 3)
		incTest(1, 0)
		incTest(3, 2)
		incTest(-2, 2)
    }
    return true
}