package cjnum.blas.testblas

import cjnum.blas.*
import cjnum.floats.*

class dsyr2kTestStruct {
    dsyr2kTestStruct(
        let ul!: Uplo,
        let tA!: Transpose,
        let n!: Int64,
        let k!: Int64,
        let alpha!: Float64,
        let a!: Array<Array<Float64>>,
        let b!: Array<Array<Float64>>,
        let c!: Array<Array<Float64>>,
        let beta!: Float64 = 0.0,
        let ans!: Array<Array<Float64>>
    ) {}
}

public func dsyr2kTest(nFloat64: NFloat64): Bool {
    let tests = [
		dsyr2kTestStruct(
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
			b: [
                [7.0, 8.0],
                [9.0, 10.0],
                [11.0, 12.0]
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
		dsyr2kTestStruct(
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
			b: [
                [7.0, 8.0],
                [9.0, 10.0],
                [11.0, 12.0]
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
		dsyr2kTestStruct(
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
			b: [
                [7.0, 8.0],
                [9.0, 10.0],
                [11.0, 12.0]
			],
			c: [
                [1.0, 2.0, 3.0],
                [0.0, 4.0, 5.0],
                [0.0, 0.0, 6.0]
			],
			beta: 2.0,
			ans: [
				[140.0, 250.0, 360.0],
                [0.0, 410.0, 568.0],
                [0.0, 0.0, 774.0]
			]
		),
		dsyr2kTestStruct(
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
			b: [
                [7.0, 8.0],
                [9.0, 10.0],
                [11.0, 12.0]
			],
			c: [
				[1.0, 0.0, 0.0],
                [2.0, 4.0, 0.0],
                [3.0, 5.0, 6.0]
			],
			beta: 2.0,
			ans: [
				[140.0, 0.0, 0.0],
                [250.0, 410.0, 0.0],
                [360.0, 568.0, 774.0]
			]
		),
		dsyr2kTestStruct(
			ul: Upper,
			tA: Trans,
			n: 3,
			k: 2,
			alpha: 3.0,
			a: [
				[1.0, 3.0, 5.0],
                [2.0, 4.0, 6.0]
			],
			b: [
				[7.0, 9.0, 11.0],
                [8.0, 10.0, 12.0]
			],
			c: [
                [1.0, 2.0, 3.0],
                [0.0, 4.0, 5.0],
                [0.0, 0.0, 6.0]
			],
			beta: 2.0,
			ans: [
				[140.0, 250.0, 360.0],
                [0.0, 410.0, 568.0],
                [0.0, 0.0, 774.0]
			]
		),
		dsyr2kTestStruct(
			ul: Lower,
			tA: Trans,
			n: 3,
			k: 2,
			alpha: 3.0,
			a: [
				[1.0, 3.0, 5.0],
                [2.0, 4.0, 6.0]
			],
			b: [
				[7.0, 9.0, 11.0],
                [8.0, 10.0, 12.0]
			],
			c: [
				[1.0, 0.0, 0.0],
                [2.0, 4.0, 0.0],
                [3.0, 5.0, 6.0]
			],
			beta: 2.0,
			ans: [
				[140.0, 0.0, 0.0],
                [250.0, 410.0, 0.0],
                [360.0, 568.0, 774.0]
			]
		),
		dsyr2kTestStruct(
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
			b: [
                [7.0, 8.0],
                [9.0, 10.0],
                [11.0, 12.0]
			],
			c: [
				[Float64.NaN, Float64.NaN, Float64.NaN],
				[-Float64.Inf, Float64.NaN, Float64.NaN],
				[-Float64.Inf, -Float64.Inf, Float64.NaN]
            ],
			ans: [
				[0.0, 0.0, 0.0],
				[-Float64.Inf, 0.0, 0.0],
				[-Float64.Inf, -Float64.Inf, 0.0]
			]
		),
		dsyr2kTestStruct(
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
			b: [
                [7.0, 8.0],
                [9.0, 10.0],
                [11.0, 12.0]
			],
			c: [
				[Float64.NaN, -Float64.Inf, -Float64.Inf],
                [Float64.NaN, Float64.NaN, -Float64.Inf],
                [Float64.NaN, Float64.NaN, Float64.NaN]
            ],
			ans: [
				[0.0, -Float64.Inf, -Float64.Inf],
                [0.0, 0.0, -Float64.Inf],
                [0.0, 0.0, 0.0]
			]
		),
		dsyr2kTestStruct(
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
			b: [
                [7.0, 8.0],
                [9.0, 10.0],
                [11.0, 12.0]
			],
			c: [
				[Float64.NaN, Float64.NaN, Float64.NaN],
				[-Float64.Inf, Float64.NaN, Float64.NaN],
				[-Float64.Inf, -Float64.Inf, Float64.NaN]
            ],
			ans: [
				[138.0, 246.0, 354.0],
                [-Float64.Inf, 402.0, 558.0],
                [-Float64.Inf, -Float64.Inf, 762.0]
			]
		),
		dsyr2kTestStruct(
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
			b: [
                [7.0, 8.0],
                [9.0, 10.0],
                [11.0, 12.0]
			],
			c: [
				[Float64.NaN, -Float64.Inf, -Float64.Inf],
                [Float64.NaN, Float64.NaN, -Float64.Inf],
                [Float64.NaN, Float64.NaN, Float64.NaN]
            ],
			ans: [
				[138.0, -Float64.Inf, -Float64.Inf],
				[246.0, 402.0, -Float64.Inf],
				[354.0, 558.0, 762.0]
			]
		),
		dsyr2kTestStruct(
			ul: Upper,
			tA: Trans,
			n: 3,
			k: 2,
			alpha: 3.0,
			a: [
				[1.0, 3.0, 5.0],
                [2.0, 4.0, 6.0]
            ],
			b: [
				[7.0, 9.0, 11.0],
                [8.0, 10.0, 12.0]
			],
			c: [
				[Float64.NaN, Float64.NaN, Float64.NaN],
				[-Float64.Inf, Float64.NaN, Float64.NaN],
				[-Float64.Inf, -Float64.Inf, Float64.NaN]
            ],
			ans: [
				[138.0, 246.0, 354.0],
                [-Float64.Inf, 402.0, 558.0],
                [-Float64.Inf, -Float64.Inf, 762.0]
			]
		),
		dsyr2kTestStruct(
			ul: Lower,
			tA: Trans,
			n: 3,
			k: 2,
			alpha: 3.0,
			a: [
				[1.0, 3.0, 5.0],
                [2.0, 4.0, 6.0]
            ],
			b: [
				[7.0, 9.0, 11.0],
                [8.0, 10.0, 12.0]
			],
			c: [
				[Float64.NaN, -Float64.Inf, -Float64.Inf],
                [Float64.NaN, Float64.NaN, -Float64.Inf],
                [Float64.NaN, Float64.NaN, Float64.NaN]
            ],
			ans: [
				[138.0, -Float64.Inf, -Float64.Inf],
				[246.0, 402.0, -Float64.Inf],
				[354.0, 558.0, 762.0]
			]
		)
    ]
	for (i in 0..tests.size) {
		let test = tests[i]
		let aFlat = flatten(test.a)
		let bFlat = flatten(test.b)
        let cFlat = flatten(test.c)
		let ansFlat = flatten(test.ans)
		nFloat64.dsyr2k(test.ul, test.tA, test.n, test.k, test.alpha, aFlat, test.a[0].size, bFlat, test.b[0].size, test.beta, cFlat, test.c[0].size)
		if (!fEqualApprox(ansFlat, cFlat, 1e-14)) {
			throw Exception("Case ${i}. Want ${ansFlat}, got ${cFlat}.")
		}
	}
    return true
}