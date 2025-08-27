package cjnum.blas.testblas

import cjnum.blas.*

class dtbmvTestStruct {
    dtbmvTestStruct(
        let ul!: Uplo,
        let tA!: Transpose,
        let d!: Diag,
        let n!: Int64,
        let k!: Int64,
        let a!: Array<Array<Float64>>,
        let x!: Array<Float64>,
        let ans!: Array<Float64>
    ) {}
}

public func dtbmvTest(nFloat64: NFloat64): Bool {
    let tests = [
        dtbmvTestStruct(
            ul: Upper,
            tA: NoTrans,
            d:  UnitDiag,
            n:  3,
            k:  1,
            a: [
                [1.0, 2.0, 0.0],
                [0.0, 1.0, 4.0],
                [0.0, 0.0, 1.0]
            ],
            x:  [2.0, 3.0, 4.0],
            ans: [8.0, 19.0, 4.0]
        ),
        dtbmvTestStruct(
            ul: Upper,
            tA: NoTrans,
            d:  NonUnitDiag,
            n:  5,
            k:  1,
            a: [
                [1.0, 3.0, 0.0, 0.0, 0.0],
                [0.0, 6.0, 7.0, 0.0, 0.0],
                [0.0, 0.0, 2.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 12.0, 3.0],
                [0.0, 0.0, 0.0, 0.0, -1.0]
            ],
            x:  [1.0, 2.0, 3.0, 4.0, 5.0],
            ans: [7.0, 33.0, 10.0, 63.0, -5.0]
        ),
        dtbmvTestStruct(
            ul: Lower,
            tA: NoTrans,
            d:  NonUnitDiag,
            n:  5,
            k:  1,
            a: [
                [7.0, 0.0, 0.0, 0.0, 0.0],
                [3.0, 6.0, 0.0, 0.0, 0.0],
                [0.0, 7.0, 2.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 12.0, 0.0],
                [0.0, 0.0, 0.0, 3.0, -1.0]
            ],
            x:  [1.0, 2.0, 3.0, 4.0, 5.0],
            ans: [7.0, 15.0, 20.0, 51.0, 7.0]
        ),
        dtbmvTestStruct(
            ul: Upper,
            tA: Trans,
            d:  NonUnitDiag,
            n:  5,
            k:  2,
            a: [
                [7.0, 3.0, 9.0, 0.0, 0.0],
                [0.0, 6.0, 7.0, 10.0, 0.0],
                [0.0, 0.0, 2.0, 1.0, 11.0],
                [0.0, 0.0, 0.0, 12.0, 3.0],
                [0.0, 0.0, 0.0, 0.0, -1.0]
            ],
            x:  [1.0, 2.0, 3.0, 4.0, 5.0],
            ans: [7.0, 15.0, 29.0, 71.0, 40.0]
        ),
        dtbmvTestStruct(
            ul: Lower,
            tA: Trans,
            d:  NonUnitDiag,
            n:  5,
            k: 2,
            a: [
                [7.0, 0.0, 0.0, 0.0, 0.0],
                [3.0, 6.0, 0.0, 0.0, 0.0],
                [9.0, 7.0, 2.0, 0.0, 0.0],
                [0.0, 10.0, 1.0, 12.0, 0.0],
                [0.0, 0.0, 11.0, 3.0, -1.0]
            ],
            x:  [1.0, 2.0, 3.0, 4.0, 5.0],
            ans: [40.0, 73.0, 65.0, 63.0, -5.0]
        )
    ]

    for (i in 0..tests.size) {
        let extra = 0
        let aFlat = if (tests[i].ul == Upper) {
            flattenBanded(tests[i].a, tests[i].k, 0)
        } else {
            flattenBanded(tests[i].a, 0, tests[i].k)
        }

        let incTest = { incX: Int64, extra: Int64 => 
            let xnew = makeIncremented(tests[i].x, incX, extra)
            let ans = makeIncremented(tests[i].ans, incX, extra)
            let lda = tests[i].k + 1
            nFloat64.dtbmv(tests[i].ul, tests[i].tA, tests[i].d, tests[i].n, tests[i].k, aFlat, lda, xnew, incX)
            if (!dSliceTolEqual(ans, xnew)) {
                throw Exception("Case ${i}, Inc ${incX}: Want ${ans}, got ${xnew}")
            }
        }

        incTest(1, extra)
        incTest(3, extra)
        incTest(-2, extra)
    }

    return true
}