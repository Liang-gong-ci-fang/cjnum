package cjnum.blas.testblas

import cjnum.blas.*
import cjnum.floats.*

class dtpsvTestStruct {
    dtpsvTestStruct(
        let n!: Int64,
        let a!: Array<Array<Float64>>,
        let ul!: Uplo,
        let tA!: Transpose,
        let d!: Diag,
        let x!: Array<Float64>,
        let ans!: Array<Float64>
    ) {}
}

public func dtpsvTest(nFloat64: NFloat64): Bool {
    let tests = [
        dtpsvTestStruct(
            n: 3,
            a: [
                [1.0, 2.0, 3.0],
                [0.0, 8.0, 15.0],
                [0.0, 0.0, 8.0]
            ],
            ul: Upper,
            tA: NoTrans,
            d: NonUnitDiag,
            x: [5.0, 6.0, 7.0],
            ans: [4.15625, -0.890625, 0.875]
        ),
        dtpsvTestStruct(
            n: 3,
            a: [
                [1.0, 2.0, 3.0],
                [0.0, 1.0, 15.0],
                [0.0, 0.0, 1.0]
            ],
            ul: Upper,
            tA: NoTrans,
            d: UnitDiag,
            x: [5.0, 6.0, 7.0],
            ans: [182.0, -99.0, 7.0]
        ),
        dtpsvTestStruct(
            n: 3,
            a: [
                [1.0, 2.0, 3.0],
                [2.0, 8.0, 0.0],
                [3.0, 15.0, 8.0]
            ],
            ul: Lower,
            tA: NoTrans,
            d: NonUnitDiag,
            x: [5.0, 6.0, 7.0],
            ans: [5.0, -0.5, -0.0625]
        ),
        dtpsvTestStruct(
            n: 3,
            a: [
                [1.0, 0.0, 0.0],
                [2.0, 8.0, 0.0],
                [3.0, 15.0, 8.0]
            ],
            ul: Lower,
            tA: NoTrans,
            d: UnitDiag,
            x: [5.0, 6.0, 7.0],
            ans: [5.0, -4.0, 52.0]
        ),
        dtpsvTestStruct(
            n: 3,
            a: [
                [1.0, 2.0, 3.0],
                [0.0, 8.0, 15.0],
                [0.0, 0.0, 8.0]
            ],
            ul: Upper,
            tA: Trans,
            d: NonUnitDiag,
            x: [5.0, 6.0, 7.0],
            ans: [5.0, -0.5, -0.0625]
        ),
        dtpsvTestStruct(
            n: 3,
            a: [
                [1.0, 2.0, 3.0],
                [0.0, 8.0, 15.0],
                [0.0, 0.0, 8.0]
            ],
            ul: Upper,
            tA: Trans,
            d: UnitDiag,
            x: [5.0, 6.0, 7.0],
            ans: [5.0, -4.0, 52.0]
        ),
        dtpsvTestStruct(
            n: 3,
            a: [
                [1.0, 0.0, 0.0],
                [2.0, 8.0, 0.0],
                [3.0, 15.0, 8.0]
            ],
            ul: Lower,
            tA: Trans,
            d: NonUnitDiag,
            x: [5.0, 6.0, 7.0],
            ans: [4.15625, -0.890625, 0.875]
        ),
        dtpsvTestStruct(
            n: 3,
            a: [
                [1.0, 0.0, 0.0],
                [2.0, 1.0, 0.0],
                [3.0, 15.0, 1.0]
            ],
            ul: Lower,
            tA: Trans,
            d: UnitDiag,
            x: [5.0, 6.0, 7.0],
            ans: [182.0, -99.0, 7.0]
        )
    ]

    for (i in 0..tests.size) {
        let incTest = { incX: Int64, extra: Int64 =>
            let aFlat = flattenTriangular(tests[i].a, tests[i].ul)
            let x = makeIncremented(tests[i].x, incX, extra)
            nFloat64.dtpsv(tests[i].ul, tests[i].tA, tests[i].d, tests[i].n, aFlat, x, incX)
            let ans = makeIncremented(tests[i].ans, incX, extra)
            if (!fEqualApprox(ans, x, 1e-14)) {
                throw Exception("Case ${i}, incX = ${incX}: Want ${ans}, got ${x}.")
            }
        }

		incTest(1, 0)
		incTest(-2, 0)
		incTest(3, 0)
		incTest(-3, 8)
		incTest(4, 2)
    }
    return true
}