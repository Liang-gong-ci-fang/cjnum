package cjnum.blas.testblas

import std.math.*
import std.random.*
import cjnum.blas.*
import cjnum.floats.*

class dsyrkTestStruct {
    dsyrkTestStruct(
        let ul!: Uplo,
        let tA!: Transpose,
        let n!: Int64,
        let k!: Int64,
        let alpha!: Float64,
        let a!: Array<Array<Float64>>,
        let c!: Array<Array<Float64>>,
        let beta!: Float64 = 0.0,
        let ans!: Array<Array<Float64>>
    ) {}
}

public func dsyrkTest(nFloat64: NFloat64): Bool {
    let tests = [
		dsyrkTestStruct(
			ul: Upper,
			tA: NoTrans,
            n: 3,
            k: 2,
            alpha: 0.0,
            a: [
				[1.0, 2.0],
				[3.0, 4.0],
                [5.0, 6.0]
            ],
            c: [
                [1.0, 2.0, 3.0],
				[0.0, 5.0, 6.0],
                [0.0, 0.0, 9.0]
			],
			beta: 2.0,
            ans: [
                [2.0, 4.0, 6.0],
                [0.0, 10.0, 12.0],
                [0.0, 0.0, 18.0]
			]
		),
		dsyrkTestStruct(
			ul: Lower,
			tA: NoTrans,
            n: 3,
            k: 2,
            alpha: 0.0,
            a: [
				[1.0, 2.0],
				[3.0, 4.0],
                [5.0, 6.0]
            ],
            c: [
                [1.0, 0.0, 0.0],
				[2.0, 3.0, 0.0],
                [4.0, 5.0, 6.0]
            ],
			beta: 2.0,
            ans: [
				[2.0, 0.0, 0.0],
				[4.0, 6.0, 0.0],
				[8.0, 10.0, 12.0]
			]
		),
		dsyrkTestStruct(
			ul: Upper,
			tA: NoTrans,
            n: 3,
            k: 2,
            alpha: 3.0,
            a: [
				[1.0, 2.0],
				[3.0, 4.0],
                [5.0, 6.0]
            ],
            c: [
				[1.0, 2.0, 3.0],
				[0.0, 4.0, 5.0],
				[0.0, 0.0, 6.0]
            ],
			beta: 2.0,
            ans: [
				[17.0, 37.0, 57.0],
				[0.0, 83.0, 127.0],
                [0.0, 0.0, 195.0]
			]
		),
		dsyrkTestStruct(
			ul: Upper,
			tA: NoTrans,
            n: 3,
            k: 2,
            alpha: 3.0,
            a: [
				[1.0, 2.0],
				[3.0, 4.0],
                [5.0, 6.0]
            ],
            c: [
				[1.0, 2.0, 3.0],
				[0.0, 4.0, 5.0],
				[0.0, 0.0, 6.0]
            ],
			beta: 0.0,
            ans: [
				[15.0, 33.0, 51.0],
				[0.0, 75.0, 117.0],
				[0.0, 0.0, 183.0]
			]
		),
		dsyrkTestStruct(
			ul: Lower,
			tA: NoTrans,
            n: 3,
            k: 2,
            alpha: 3.0,
            a: [
				[1.0, 2.0],
				[3.0, 4.0],
                [5.0, 6.0]
            ],
            c: [
                [1.0, 0.0, 0.0],
				[2.0, 4.0, 0.0],
                [3.0, 5.0, 6.0]
            ],
			beta: 2.0,
            ans: [
				[17.0, 0.0, 0.0],
                [37.0, 83.0, 0.0],
                [57.0, 127.0, 195.0]
			]
		),
		dsyrkTestStruct(
			ul: Lower,
			tA: NoTrans,
            n: 3,
            k: 2,
            alpha: 3.0,
            a: [
				[1.0, 2.0],
				[3.0, 4.0],
                [5.0, 6.0]
            ],
            c: [
                [1.0, 0.0, 0.0],
				[2.0, 4.0, 0.0],
                [3.0, 5.0, 6.0]
            ],
			beta: 0.0,
            ans: [
				[15.0, 0.0, 0.0],
				[33.0, 75.0, 0.0],
                [51.0, 117.0, 183.0]
			]
		),
		dsyrkTestStruct(
			ul: Upper,
			tA: Trans,
            n: 3,
            k: 2,
            alpha: 3.0,
            a: [
				[1.0, 3.0, 5.0],
				[2.0, 4.0, 6.0]
            ],
            c: [
				[1.0, 2.0, 3.0],
				[0.0, 4.0, 5.0],
				[0.0, 0.0, 6.0]
            ],
			beta: 2.0,
            ans: [
				[17.0, 37.0, 57.0],
				[0.0, 83.0, 127.0],
				[0.0, 0.0, 195.0]
			]
		),
		dsyrkTestStruct(
			ul: Lower,
			tA: Trans,
            n: 3,
            k: 2,
            alpha: 3.0,
            a: [
				[1.0, 3.0, 5.0],
				[2.0, 4.0, 6.0]
            ],
            c: [
				[1.0, 0.0, 0.0],
				[2.0, 4.0, 0.0],
				[3.0, 5.0, 6.0]
            ],
			beta: 2.0,
            ans: [
				[17.0, 0.0, 0.0],
				[37.0, 83.0, 0.0],
				[57.0, 127.0, 195.0]
			]
		),
		dsyrkTestStruct(
			ul: Upper,
			tA: Trans,
            n: 3,
            k: 2,
            alpha: 3.0,
            a: [
				[1.0, 3.0, 5.0],
				[2.0, 4.0, 6.0]
            ],
            c: [
				[1.0, 2.0, 3.0],
				[0.0, 4.0, 5.0],
				[0.0, 0.0, 6.0]
            ],
			beta: 0.0,
            ans: [
				[15.0, 33.0, 51.0],
				[0.0, 75.0, 117.0],
				[0.0, 0.0, 183.0]
			]
		),
		dsyrkTestStruct(
			ul: Lower,
			tA: Trans,
            n: 3,
            k: 2,
            alpha: 3.0,
            a: [
				[1.0, 3.0, 5.0],
				[2.0, 4.0, 6.0]
            ],
            c: [
				[1.0, 0.0, 0.0],
				[2.0, 4.0, 0.0],
				[3.0, 5.0, 6.0]
            ],
			beta: 0.0,
            ans: [
				[15.0, 0.0, 0.0],
				[33.0, 75.0, 0.0],
				[51.0, 117.0, 183.0]
			]
		),
		dsyrkTestStruct(
			ul: Lower,
			tA: Trans,
            n: 3,
            k: 2,
            alpha: 0.0,
            a: [
				[1.0, 3.0, 5.0],
				[2.0, 4.0, 6.0]
            ],
            c: [
				[1.0, 0.0, 0.0],
				[2.0, 4.0, 0.0],
				[3.0, 5.0, 6.0]
            ],
			beta: 0.0,
            ans: [
				[0.0, 0.0, 0.0],
				[0.0, 0.0, 0.0],
				[0.0, 0.0, 0.0]
			]
		),
		dsyrkTestStruct(
			ul: Upper,
			tA: Trans,
            n: 3,
            k: 2,
            alpha: 0.0,
            a: [
				[1.0, 3.0, 5.0],
				[2.0, 4.0, 6.0]
            ],
            c: [
				[1.0, 2.0, 3.0],
				[0.0, 4.0, 5.0],
				[0.0, 0.0, 6.0]
            ],
			beta: 0.0,
            ans: [
				[0.0, 0.0, 0.0],
				[0.0, 0.0, 0.0],
				[0.0, 0.0, 0.0]
			]
		)
    ]
    for (i in 0..tests.size) {
        let test = tests[i]
        let aFlat = flatten(test.a)
        let cFlat = flatten(test.c)
        let ansFlat = flatten(test.ans)
        nFloat64.dsyrk(test.ul, test.tA, test.n, test.k, test.alpha, aFlat, test.a[0].size, test.beta, cFlat, test.c[0].size)
        if (!fEqualApprox(ansFlat, cFlat, 1e-14)) {
            throw Exception("Case ${i}. Want ${ansFlat}, got ${cFlat}.")
        }
    }
    return true
}