package cjnum.blas.testblas

import cjnum.blas.*
import cjnum.floats.*

class dsbmvTestStruct {
    dsbmvTestStruct(
        let ul!: Uplo,
        let n!: Int64,
        let k!: Int64,
        let alpha!: Float64,
        let beta!: Float64,
        let a!: Array<Array<Float64>>,
        let x!: Array<Float64>,
        let y!: Array<Float64>,
        let ans!: Array<Float64>
    ) {}
}

public func dsbmvTest(nFloat64: NFloat64): Bool {
    let tests = [
        dsbmvTestStruct(
            ul: Upper,
            n: 4, 
            k: 2,
            alpha: 2.0,
            beta: 3.0,
            a: [
                [7.0, 8.0, 2.0, 0.0],
                [0.0, 8.0, 2.0, -3.0],
                [0.0, 0.0, 3.0, 6.0],
                [0.0, 0.0, 0.0, 9.0]
            ],
            x: [1.0, 2.0, 3.0, 4.0],
            y: [-1.0, -2.0, -3.0, -4.0],
            ans: [55.0, 30.0, 69.0, 84.0]
        ),
        dsbmvTestStruct(
            ul: Lower,
            n: 4, 
            k: 2,
            alpha: 2.0,
            beta: 3.0,
            a: [
                [7.0, 0.0, 0.0, 0.0],
                [8.0, 8.0, 0.0, 0.0],
                [2.0, 2.0, 3.0, 0.0],
                [0.0, -3.0, 6.0, 9.0]
            ],
            x: [1.0, 2.0, 3.0, 4.0],
            y: [-1.0, -2.0, -3.0, -4.0],
            ans: [55.0, 30.0, 69.0, 84.0]
        )
    ]

    for (i in 0..tests.size) {
        let extra = 0
        let aFlat = if (tests[i].ul == Upper) {
            flattenBanded(tests[i].a, tests[i].k, 0)
        } else {
            flattenBanded(tests[i].a, 0, tests[i].k)
        }
        let incTest = { incX: Int64, incY: Int64, extra: Int64 =>
            let xnew = makeIncremented(tests[i].x, incX, extra)
            let ynew = makeIncremented(tests[i].y, incY, extra)
            let ans = makeIncremented(tests[i].ans, incY, extra)
            nFloat64.dsbmv(tests[i].ul, tests[i].n, tests[i].k, tests[i].alpha, aFlat, tests[i].k + 1, xnew, incX, tests[i].beta, ynew, incY)
            if (!dSliceTolEqual(ans, ynew)) {
                throw Exception("Case ${i} incX = ${incX} incY = ${incY}: Want ${ans}, got ${ynew}")
            }
        }

		incTest(1, 1, extra)
		incTest(1, 3, extra)
		incTest(1, -3, extra)
		incTest(2, 3, extra)
		incTest(2, -3, extra)
		incTest(3, 2, extra)
		incTest(-3, 2, extra)
    }
    return true
}