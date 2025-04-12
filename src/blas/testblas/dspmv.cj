package cjnum.blas.testblas

import cjnum.blas.*
import cjnum.floats.*

class dspmvTestStruct {
    dspmvTestStruct(
        let ul!: Uplo,
        let n!: Int64,
        let a!: Array<Array<Float64>>,
        let x!: Array<Float64>,
        let y!: Array<Float64>,
        let alpha!: Float64,
        let beta!: Float64,
        let ans!: Array<Float64>
    ) {}
}

public func dspmvTest(nFloat64: NFloat64): Bool {
    let tests = [
        dspmvTestStruct(
            ul: Upper,
            n: 3,
            a: [
                [5.0, 6.0, 7.0],
                [0.0, 8.0, 10.0],
                [0.0, 0.0, 13.0]
            ],
            x: [3.0, 4.0, 5.0],
            y: [6.0, 7.0, 8.0],
            alpha: 2.1,
            beta: -3.0,
            ans: [137.4, 189.0, 240.6]
        ),
        dspmvTestStruct(
            ul: Lower,
            n: 3,
            a: [
                [5.0, 0.0, 0.0],
                [6.0, 8.0, 0.0],
                [7.0, 10.0, 13.0]
            ],
            x: [3.0, 4.0, 5.0],
            y: [6.0, 7.0, 8.0],
            alpha: 2.1,
            beta: -3.0,
            ans: [137.4, 189.0, 240.6]
        )
    ]

    for (i in 0..tests.size) {
        let incTest = { incX: Int64, incY: Int64, extra: Int64 =>
            let x = makeIncremented(tests[i].x, incX, extra)
            let y = makeIncremented(tests[i].y, incY, extra)
            let aFlat = flattenTriangular(tests[i].a, tests[i].ul)
            let ans = makeIncremented(tests[i].ans, incY, extra)

            nFloat64.dspmv(tests[i].ul, tests[i].n, tests[i].alpha, aFlat, x, incX, tests[i].beta, y, incY)
            if (!fEqualApprox(ans, y, 1e-14)) {
                throw Exception("Case ${i}, incX=${incX}, incY=${incY}: Want ${ans}, got ${y}.")
            }
        }
		incTest(1, 1, 0)
		incTest(2, 3, 0)
		incTest(3, 2, 0)
		incTest(-3, 2, 0)
		incTest(-2, 4, 0)
		incTest(2, -1, 0)
		incTest(-3, -4, 3)
    }
    return true
}