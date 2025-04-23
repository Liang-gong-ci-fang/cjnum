package cjnum.blas.testblas

import cjnum.blas.*
import cjnum.floats.*

class dspr2TestStruct {
    dspr2TestStruct(
        let n!: Int64,
        let a!: Array<Array<Float64>>,
        let ul!: Uplo,
        let x!: Array<Float64>,
        let y!: Array<Float64>,
        let alpha!: Float64,
        let ans!: Array<Array<Float64>>
    ) {}
}

public func dspr2Test(nFloat64: NFloat64): Bool {
    let tests = [
        dspr2TestStruct(
            n: 3, 
            a: [
                [7.0, 2.0, 4.0],
                [0.0, 3.0, 5.0],
                [0.0, 0.0, 6.0]
            ],
            x: [2.0, 3.0, 4.0],
            y: [5.0, 6.0, 7.0],
            alpha: 2.0,
            ul: Upper,
            ans: [
                [47.0, 56.0, 72.0],
                [0.0, 75.0, 95.0],
                [0.0, 0.0, 118.0]
            ]
        ),
        dspr2TestStruct(
            n: 3,
            a: [
                [7.0, 0.0, 0.0],
                [2.0, 3.0, 0.0],
                [4.0, 5.0, 6.0]
            ],
            x: [2.0, 3.0, 4.0],
            y: [5.0, 6.0, 7.0],
            alpha: 2.0,
            ul: Lower,
            ans: [
                [47.0, 0.0, 0.0],
                [56.0, 75.0, 0.0],
                [72.0, 95.0, 118.0]
            ]
        )
    ]
    
    for (i in 0..tests.size) {
        let incTest = { incX: Int64, incY: Int64, extra: Int64 =>
            let aFlat = flattenTriangular(tests[i].a, tests[i].ul)
            let x = makeIncremented(tests[i].x, incX, extra)
            let y = makeIncremented(tests[i].y, incY, extra)
            nFloat64.dspr2(tests[i].ul, tests[i].n, tests[i].alpha, x, incX, y, incY, aFlat)
            let ansFlat = flattenTriangular(tests[i].ans, tests[i].ul)
            if (!fEqualApprox(aFlat, ansFlat, 1e-14)) {
                throw Exception("Case${i}, incX = ${incX}, incY = ${incY}. Want ${ansFlat}, got ${aFlat}.")
            }
        }
		incTest(1, 1, 0)
		incTest(-2, 1, 0)
		incTest(-2, 3, 0)
		incTest(2, -3, 0)
		incTest(3, -2, 0)
		incTest(-3, -4, 0)
    }
    return true
}