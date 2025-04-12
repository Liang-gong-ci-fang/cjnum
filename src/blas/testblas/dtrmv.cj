package cjnum.blas.testblas

import cjnum.blas.*
import cjnum.floats.*

class dtrmvTestStruct {
    dtrmvTestStruct(
        let n!: Int64,
        let a!: Array<Array<Float64>>,
        let x!: Array<Float64>,
        let d!: Diag,
        let ul!: Uplo,
        let tA!: Transpose,
        let ans!: Array<Float64>
    ) {}
}

public func dtrmvTest(nFloat64: NFloat64): Bool {
    let tests = [
        dtrmvTestStruct(
            n: 1,
            a: [[5.0]],
            x: [2.0],
            d: NonUnitDiag,
            ul: Upper,
            tA: NoTrans,
            ans: [10.0]
        ),
        dtrmvTestStruct(
            n: 1,
            a: [[5.0]],
            x: [2.0],
            d: UnitDiag,
            ul: Upper,
            tA: NoTrans,
            ans: [2.0]
        ),
        dtrmvTestStruct(
            n: 3,
            a: [[5.0, 6.0, 7.0], [0.0, 9.0, 10.0], [0.0, 0.0, 13.0]],
            x: [3.0, 4.0, 5.0],
            d: NonUnitDiag,
            ul: Upper,
            tA: NoTrans,
            ans: [74.0, 86.0, 65.0]
        ),
        dtrmvTestStruct(
            n: 3,
            a: [[5.0, 6.0, 7.0], [0.0, 9.0, 10.0], [0.0, 0.0, 13.0]],
            x: [3.0, 4.0, 5.0],
            d: UnitDiag,
            ul: Upper,
            tA: NoTrans,
            ans: [62.0, 54.0, 5.0]
        ),
        dtrmvTestStruct(
            n: 3,
            a: [[5.0, 0.0, 0.0], [6.0, 9.0, 0.0], [7.0, 10.0, 13.0]],
            x: [3.0, 4.0, 5.0],
            d: NonUnitDiag,
            ul: Lower,
            tA: NoTrans,
            ans: [15.0, 54.0, 126.0]
        ),
        dtrmvTestStruct(
            n: 3,
            a: [[1.0, 0.0, 0.0], [6.0, 1.0, 0.0], [7.0, 10.0, 1.0]],
            x: [3.0, 4.0, 5.0],
            d: UnitDiag,
            ul: Lower,
            tA: NoTrans,
            ans: [3.0, 22.0, 66.0]
        ),
        dtrmvTestStruct(
            n: 3,
            a: [[5.0, 6.0, 7.0], [0.0, 9.0, 10.0], [0.0, 0.0, 13.0]],
            x: [3.0, 4.0, 5.0],
            d: NonUnitDiag,
            ul: Upper,
            tA: Trans,
            ans: [15.0, 54.0, 126.0]
        ),
        dtrmvTestStruct(
            n: 3,
            a: [[1.0, 6.0, 7.0], [0.0, 1.0, 10.0], [0.0, 0.0, 1.0]],
            x: [3.0, 4.0, 5.0],
            d: UnitDiag,
            ul: Upper,
            tA: Trans,
            ans: [3.0, 22.0, 66.0]
        ),
        dtrmvTestStruct(
            n: 3,
            a: [[5.0, 0.0, 0.0], [6.0, 9.0, 0.0], [7.0, 10.0, 13.0]],
            x: [3.0, 4.0, 5.0],
            d: NonUnitDiag,
            ul: Lower,
            tA: Trans,
            ans: [74.0, 86.0, 65.0]
        )
    ]

    for (i in 0..tests.size) {
        let incTest = { incX: Int64, extra: Int64 => 
            let aFlat = flatten(tests[i].a)
            let x = makeIncremented(tests[i].x, incX, extra)
            nFloat64.dtrmv(tests[i].ul, tests[i].tA, tests[i].d, tests[i].n, aFlat, tests[i].n, x, incX)
            let ans = makeIncremented(tests[i].ans, incX, extra)
            if (!fEqualApprox(x, ans, 1e-14)) {
                throw Exception("Case ${i}, idx ${incX}  Want ${ans}, got ${x}.")
            }
        }

        incTest(1, 3)
        incTest(-3, 3)
        incTest(4, 3)
    }

    return true
}
