package cjnum.blas.testblas

import cjnum.blas.*

class dgbmvTestStruct {
    dgbmvTestStruct(
        let tA !: Transpose,
        let m!: Int64,
        let n!: Int64,
        let kL !: Int64,
        let kU!: Int64,
        let alpha!: Float64,
        let a!: Array<Array<Float64>>,
        let lda!: Int64,
        let x!: Array<Float64>,
        let beta!: Float64,
        let y!: Array<Float64>,
        let ans!: Array<Float64>
    ) {}
}

public func dgbmvTest(nFloat64: NFloat64): Bool {
    let tests = [
        dgbmvTestStruct(
            tA:    NoTrans,
            m:     9,
            n:     6,
            lda:   4,
            kL:    2,
            kU:    1,
            alpha: 3.0,
            beta:  2.0,
            a: [
                [5.0, 3.0, 0.0, 0.0, 0.0, 0.0],
                [-1.0, 2.0, 9.0, 0.0, 0.0, 0.0],
                [4.0, 8.0, 3.0, 6.0, 0.0, 0.0],
                [0.0, -1.0, 8.0, 2.0, 1.0, 0.0],
                [0.0, 0.0, 9.0, 9.0, 9.0, 5.0],
                [0.0, 0.0, 0.0, 2.0, -3.0, 2.0],
                [0.0, 0.0, 0.0, 0.0, 1.0, 5.0],
                [0.0, 0.0, 0.0, 0.0, 0.0, 6.0]
            ],
            x: [1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
            y: [-1.0, -2.0, -3.0, -4.0, -5.0, -6.0, -7.0, -8.0, -9.0],
            ans: [31.0, 86.0, 153.0, 97.0, 404.0, 3.0, 91.0, 92.0, -18.0]
        ),
        dgbmvTestStruct(
            tA:    Trans,
            m:     9,
            n:     6,
            lda:   4,
            kL:    2,
            kU:    1,
            alpha: 3.0,
            beta:  2.0,
            a: [
                [5.0, 3.0, 0.0, 0.0, 0.0, 0.0],
                [-1.0, 2.0, 9.0, 0.0, 0.0, 0.0],
                [4.0, 8.0, 3.0, 6.0, 0.0, 0.0],
                [0.0, -1.0, 8.0, 2.0, 1.0, 0.0],
                [0.0, 0.0, 9.0, 9.0, 9.0, 5.0],
                [0.0, 0.0, 0.0, 2.0, -3.0, 2.0],
                [0.0, 0.0, 0.0, 0.0, 1.0, 5.0],
                [0.0, 0.0, 0.0, 0.0, 0.0, 6.0]
            ],
            x: [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0],
            y: [-1.0, -2.0, -3.0, -4.0, -5.0, -6.0],
            ans: [43.0, 77.0, 306.0, 241.0, 104.0, 348.0]
        ),
        dgbmvTestStruct(
            tA:    NoTrans,
            m:     6,
            n:     3,
            lda:   1,
            kL:    0,
            kU:    0,
            alpha: 2.0,
            beta:  1.0,
            a: [
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0]
            ],
            x: [1.0, 2.0, 3.0],
            y: [-1.0, -2.0, -3.0, -4.0, -5.0, -6.0],
            ans: [1.0, 2.0, 3.0, -4.0, -5.0, -6.0]
        ),
        dgbmvTestStruct(
            tA:    Trans,
            m:     6,
            n:     3,
            lda:   1,
            kL:    0,
            kU:    0,
            alpha: 2.0,
            beta:  1.0,
            a: [
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
                [0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0]
            ],
            x: [1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
            y: [-1.0, -2.0, -3.0],
            ans: [1.0, 2.0, 3.0]
        ),
        dgbmvTestStruct(
            tA:    NoTrans,
            m:     3,
            n:     5,
            lda:   4,
            kL:    1,
            kU:    2,
            alpha: 2.0,
            beta:  1.0,
            a: [
                [1.0, 2.0, 3.0, 0.0, 0.0],
                [1.0, 3.0, 6.0, 9.0, 0.0],
                [0.0, 1.0, 1.0, 1.0, 1.0]
            ],
            x: [1.0, 2.0, 3.0, 4.0, 5.0],
            y: [-1.0, -2.0, -3.0],
            ans: [-1.0 + 2.0 * (1.0 + 4.0 + 9.0), -2.0 + 2.0 * (1.0 + 6.0 + 18.0 + 36.0), -3.0 + 2.0 * (2.0 + 3.0 + 4.0 + 5.0)]
        )
    ]
    for (i in 0..tests.size) {
        let extra = 3
        let aFlat = flattenBanded(tests[i].a, tests[i].kU, tests[i].kL)
        let incTest = { incX: Int64, incY: Int64, extra: Int64 =>
            let xnew = makeIncremented(tests[i].x, incX, extra)
            let ynew = makeIncremented(tests[i].y, incY, extra)
            let ans = makeIncremented(tests[i].ans, incY, extra)
            nFloat64.dgbmv(tests[i].tA, tests[i].m, tests[i].n, tests[i].kL, tests[i].kU, tests[i].alpha, aFlat, tests[i].lda, xnew, incX, tests[i].beta, ynew, incY)
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