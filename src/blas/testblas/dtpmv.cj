package cjnum.blas.testblas

import cjnum.blas.*
import cjnum.floats.*

class dtpmvTestStruct {
    dtpmvTestStruct(
        let n!: Int64,
        let a!: Array<Array<Float64>>,
        let x!: Array<Float64>,
        let d!: Diag,
        let ul!: Uplo,
        let tA!: Transpose,
        let ans!: Array<Float64>
    ) {}
}

public func dtpmvTest(nFLoat64: NFloat64): Bool {
    let tests = [
        dtpmvTestStruct(
            n: 3,
            a: [
                [5.0, 6.0, 7.0],
                [0.0, 9.0, 10.0],
                [0.0, 0.0, 13.0]
            ],
            x: [3.0, 4.0, 5.0],
            d: NonUnitDiag,
            ul: Upper,
            tA: NoTrans,
            ans: [74.0, 86.0, 65.0]
        ),
        dtpmvTestStruct(
            n: 3,
            a: [
                [5.0, 6.0, 7.0],
                [0.0, 9.0, 10.0],
                [0.0, 0.0, 13.0]
            ],
            x: [3.0, 4.0, 5.0],
            d: UnitDiag,
            ul: Upper,
            tA: NoTrans,
            ans: [62.0, 54.0, 5.0]
        ),
        dtpmvTestStruct(
            n: 3,
            a: [
                [5.0, 0.0, 0.0],
                [6.0, 9.0, 0.0],
                [7.0, 10.0, 13.0]
            ],
            x: [3.0, 4.0, 5.0],
            d: NonUnitDiag,
            ul: Lower,
            tA: NoTrans,
            ans: [15.0, 54.0, 126.0]
        ),
        dtpmvTestStruct(
            n: 3,
            a: [
                [1.0, 0.0, 0.0],
                [6.0, 1.0, 0.0],
                [7.0, 10.0, 1.0]
            ],
            x: [3.0, 4.0, 5.0],
            d: UnitDiag,
            ul: Lower,
            tA: NoTrans,
            ans: [3.0, 22.0, 66.0]
        ),
        dtpmvTestStruct(
            n: 3,
            a: [
                [5.0, 6.0, 7.0],
                [0.0, 9.0, 10.0],
                [0.0, 0.0, 13.0]
            ],
            x: [3.0, 4.0, 5.0],
            d: NonUnitDiag,
            ul: Upper,
            tA: Trans,
            ans: [15.0, 54.0, 126.0]
        ),
        dtpmvTestStruct(
            n: 3,
            a: [
                [1.0, 6.0, 7.0],
                [0.0, 1.0, 10.0],
                [0.0, 0.0, 1.0]
            ],
            x: [3.0, 4.0, 5.0],
            d: UnitDiag,
            ul: Upper,
            tA: Trans,
            ans: [3.0, 22.0, 66.0]
        ),
        dtpmvTestStruct(
            n: 3,
            a: [
                [5.0, 0.0, 0.0],
                [6.0, 9.0, 0.0],
                [7.0, 10.0, 13.0]
            ],
            x: [3.0, 4.0, 5.0],
            d: NonUnitDiag,
            ul: Lower,
            tA: Trans,
            ans: [74.0, 86.0, 65.0]
        )
    ]

    for (i in 0..tests.size) {
        let incTest = { incX: Int64, extra: Int64 =>
            let aFlat = flattenTriangular(tests[i].a, tests[i].ul)
            let x = makeIncremented(tests[i].x, incX, extra)
            nFLoat64.dtpmv(tests[i].ul, tests[i].tA, tests[i].d, tests[i].n, aFlat, x, incX)
            let ans = makeIncremented(tests[i].ans, incX, extra)
            if (!fEqualApprox(x, ans, 1e-14)) {
                throw Exception("Case ${i}, idx ${incX}: Want ${ans}, got ${x}.")
            }
        }

        incTest(1, 0)
        incTest(-3, 3)
        incTest(4, 3)
    }
    return true
}