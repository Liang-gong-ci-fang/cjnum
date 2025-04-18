package cjnum.blas.testblas

import std.math.*
import std.random.*
import cjnum.blas.*
import cjnum.floats.*

class dsymmTestStruct {
    dsymmTestStruct(
        let m!: Int64,
        let n!: Int64,
        let side!: Side,
        let ul!: Uplo,
        let a!: Array<Array<Float64>>,
        let b!: Array<Array<Float64>>,
        let c!: Array<Array<Float64>>,
        let alpha!: Float64,
        let beta!: Float64 = 0.0,
        let ans!: Array<Array<Float64>>
    ) {}
}

public func dsymmTest(nFloat64: NFloat64): Bool {
    let tests = [
		dsymmTestStruct(
			side: Left,
			ul: Upper,
			m: 3,
			n: 4,
			a: [
				[2.0, 3.0, 4.0],
				[0.0, 6.0, 7.0],
                [0.0, 0.0, 10.0]
            ],
            b: [
                [2.0, 3.0, 4.0, 8.0],
				[5.0, 6.0, 7.0, 15.0],
				[8.0, 9.0, 10.0, 20.0]
            ],
            c: [
                [8.0, 12.0, 2.0, 1.0],
                [9.0, 12.0, 9.0, 9.0],
				[12.0, 1.0, -1.0, 5.0]
            ],
            alpha: 2.0,
            beta: 3.0,
            ans: [
				[126.0, 156.0, 144.0, 285.0],
                [211.0, 252.0, 275.0, 535.0],
                [282.0, 291.0, 327.0, 689.0]
			]
		),
		dsymmTestStruct(
            side: Left,
            ul: Upper,
            m: 3,
            n: 4,
			a: [
                [0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0]
			],
			b: [
				[0.0, 0.0, 0.0, 0.0],
				[0.0, 0.0, 0.0, 0.0],
				[0.0, 0.0, 0.0, 0.0]
            ],
            c: [
				[8.0, 12.0, 2.0, 1.0],
				[9.0, 12.0, 9.0, 9.0],
				[12.0, 1.0, -1.0, 5.0]
			],
			alpha: 0.0,
			beta: 2.0,
            ans: [
                [16.0, 24.0, 4.0, 2.0],
                [18.0, 24.0, 18.0, 18.0],
                [24.0, 2.0, -2.0, 10.0]
			]
		),
		dsymmTestStruct(
            side: Left,
            ul: Upper,
            m: 4,
            n: 3,
			a: [
				[2.0, 3.0, 4.0, 8.0],
				[0.0, 6.0, 7.0, 9.0],
                [0.0, 0.0, 10.0, 10.0],
                [0.0, 0.0, 0.0, 11.0]
			],
			b: [
				[2.0, 3.0, 4.0],
                [5.0, 6.0, 7.0],
                [8.0, 9.0, 10.0],
                [2.0, 1.0, 1.0]
			],
			c: [
				[8.0, 12.0, 2.0],
				[9.0, 12.0, 9.0],
				[12.0, 1.0, -1.0],
				[1.0, 9.0, 5.0]
            ],
            alpha: 2.0,
            beta: 3.0,
            ans: [
				[158.0, 172.0, 160.0],
				[247.0, 270.0, 293.0],
                [322.0, 311.0, 347.0],
                [329.0, 385.0, 427.0]
			]
		),
		dsymmTestStruct(
            side: Left,
            ul: Lower,
            m: 3,
            n: 4,
			a: [
                [2.0, 0.0, 0.0],
                [3.0, 6.0, 0.0],
                [4.0, 7.0, 10.0]
            ],
			b: [
                [2.0, 3.0, 4.0, 8.0],
                [5.0, 6.0, 7.0, 15.0],
                [8.0, 9.0, 10.0, 20.0]
            ],
			c: [
                [8.0, 12.0, 2.0, 1.0],
                [9.0, 12.0, 9.0, 9.0],
                [12.0, 1.0, -1.0, 5.0]
            ],
			alpha: 2.0,
            beta: 3.0,
            ans: [
                [126.0, 156.0, 144.0, 285.0],
                [211.0, 252.0, 275.0, 535.0],
				[282.0, 291.0, 327.0, 689.0]
			]
		),
		dsymmTestStruct(
			side: Left,
			ul: Lower,
			m: 4,
			n: 3,
			a: [
				[2.0, 0.0, 0.0, 0.0],
				[3.0, 6.0, 0.0, 0.0],
                [4.0, 7.0, 10.0, 0.0],
                [8.0, 9.0, 10.0, 11.0]
			],
			b: [
                [2.0, 3.0, 4.0],
                [5.0, 6.0, 7.0],
                [8.0, 9.0, 10.0],
				[2.0, 1.0, 1.0]
			],
			c: [
				[8.0, 12.0, 2.0],
                [9.0, 12.0, 9.0],
                [12.0, 1.0, -1.0],
                [1.0, 9.0, 5.0]
			],
			alpha: 2.0,
            beta: 3.0,
            ans: [
				[158.0, 172.0, 160.0],
				[247.0, 270.0, 293.0],
				[322.0, 311.0, 347.0],
				[329.0, 385.0, 427.0]
			]
		),
		dsymmTestStruct(
            side: Right,
            ul: Upper,
            m: 3,
            n: 4,
			a: [
				[2.0, 0.0, 0.0, 0.0],
				[3.0, 6.0, 0.0, 0.0],
                [4.0, 7.0, 10.0, 0.0],
                [3.0, 4.0, 5.0, 6.0]
			],
			b: [
                [2.0, 3.0, 4.0, 9.0],
                [5.0, 6.0, 7.0, -3.0],
                [8.0, 9.0, 10.0, -2.0]
            ],
			c: [
                [8.0, 12.0, 2.0, 10.0],
                [9.0, 12.0, 9.0, 10.0],
                [12.0, 1.0, -1.0, 10.0]
            ],
			alpha: 2.0,
            beta: 3.0,
            ans: [
                [32.0, 72.0, 86.0, 138.0],
                [47.0, 108.0, 167.0, -6.0],
				[68.0, 111.0, 197.0, 6.0]
			]
		),
		dsymmTestStruct(
            side: Right,
            ul: Upper,
            m: 4,
            n: 3,
			a: [
                [2.0, 0.0, 0.0],
                [3.0, 6.0, 0.0],
                [4.0, 7.0, 10.0]
			],
			b: [
                [2.0, 3.0, 4.0],
                [5.0, 6.0, 7.0],
                [8.0, 9.0, 10.0],
				[2.0, 1.0, 1.0]
            ],
			c: [
				[8.0, 12.0, 2.0],
                [9.0, 12.0, 9.0],
                [12.0, 1.0, -1.0],
                [1.0, 9.0, 5.0]
			],
			alpha: 2.0,
            beta: 3.0,
            ans: [
                [32.0, 72.0, 86.0],
                [47.0, 108.0, 167.0],
                [68.0, 111.0, 197.0],
				[11.0, 39.0, 35.0]
			]
		),
		dsymmTestStruct(
            side: Right,
            ul: Lower,
            m: 3,
            n: 4,
			a: [
                [2.0, 0.0, 0.0, 0.0],
                [3.0, 6.0, 0.0, 0.0],
                [4.0, 7.0, 10.0, 0.0],
				[3.0, 4.0, 5.0, 6.0]
			],
			b: [
                [2.0, 3.0, 4.0, 2.0],
                [5.0, 6.0, 7.0, 1.0],
                [8.0, 9.0, 10.0, 1.0]
            ],
			c: [
                [8.0, 12.0, 2.0, 1.0],
                [9.0, 12.0, 9.0, 9.0],
                [12.0, 1.0, -1.0, 5.0]
            ],
			alpha: 2.0,
            beta: 3.0,
            ans: [
                [94.0, 156.0, 164.0, 103.0],
                [145.0, 244.0, 301.0, 187.0],
				[208.0, 307.0, 397.0, 247.0]
			]
		),
		dsymmTestStruct(
			side: Right,
			ul: Lower,
			m: 4,
			n: 3,
			a: [
				[2.0, 0.0, 0.0],
				[3.0, 6.0, 0.0],
                [4.0, 7.0, 10.0]
			],
			b: [
                [2.0, 3.0, 4.0],
                [5.0, 6.0, 7.0],
                [8.0, 9.0, 10.0],
				[2.0, 1.0, 1.0]
			],
			c: [
				[8.0, 12.0, 2.0],
                [9.0, 12.0, 9.0],
                [12.0, 1.0, -1.0],
                [1.0, 9.0, 5.0]
			],
			alpha: 2.0,
            beta: 3.0,
            ans: [
				[82.0, 140.0, 144.0],
				[139.0, 236.0, 291.0],
                [202.0, 299.0, 387.0],
                [25.0, 65.0, 65.0]
			]
		),
		dsymmTestStruct(
            side: Left,
            ul: Upper,
            m: 3,
            n: 4,
			a: [
				[2.0, 3.0, 4.0],
                [0.0, 6.0, 7.0],
                [0.0, 0.0, 10.0]
			],
			b: [
                [2.0, 3.0, 4.0, 8.0],
                [5.0, 6.0, 7.0, 15.0],
                [8.0, 9.0, 10.0, 20.0]
            ],
			c: [
				[Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN],
				[Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN],
				[Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN]
			],
			alpha: 2.0,
            ans: [
                [102.0, 120.0, 138.0, 282.0],
                [184.0, 216.0, 248.0, 508.0],
				[246.0, 288.0, 330.0, 674.0]
			]
		),
		dsymmTestStruct(
            side: Left,
            ul: Upper,
            m: 4,
            n: 3,
			a: [
				[2.0, 3.0, 4.0, 8.0],
				[0.0, 6.0, 7.0, 9.0],
                [0.0, 0.0, 10.0, 10.0],
                [0.0, 0.0, 0.0, 11.0]
			],
			b: [
                [2.0, 3.0, 4.0],
                [5.0, 6.0, 7.0],
                [8.0, 9.0, 10.0],
				[2.0, 1.0, 1.0]
            ],
			c: [
                [Float64.NaN, Float64.NaN, Float64.NaN],
                [Float64.NaN, Float64.NaN, Float64.NaN],
                [Float64.NaN, Float64.NaN, Float64.NaN],
				[Float64.NaN, Float64.NaN, Float64.NaN]
            ],
			alpha: 2.0,
            ans: [
                [134.0, 136.0, 154.0],
                [220.0, 234.0, 266.0],
                [286.0, 308.0, 350.0],
				[326.0, 358.0, 412.0]
			]
		),
		dsymmTestStruct(
            side: Left,
            ul: Lower,
            m: 3,
            n: 4,
			a: [
				[2.0, 0.0, 0.0],
                [3.0, 6.0, 0.0],
                [4.0, 7.0, 10.0]
			],
			b: [
                [2.0, 3.0, 4.0, 8.0],
                [5.0, 6.0, 7.0, 15.0],
                [8.0, 9.0, 10.0, 20.0]
            ],
			c: [
                [Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN],
                [Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN],
                [Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN]
            ],
			alpha: 2.0,
            ans: [
                [102.0, 120.0, 138.0, 282.0],
                [184.0, 216.0, 248.0, 508.0],
                [246.0, 288.0, 330.0, 674.0]
			]
		),
		dsymmTestStruct(
            side: Left,
            ul: Lower,
            m: 4,
            n: 3,
			a: [
				[2.0, 0.0, 0.0, 0.0],
				[3.0, 6.0, 0.0, 0.0],
                [4.0, 7.0, 10.0, 0.0],
                [8.0, 9.0, 10.0, 11.0]
			],
			b: [
                [2.0, 3.0, 4.0],
                [5.0, 6.0, 7.0],
                [8.0, 9.0, 10.0],
				[2.0, 1.0, 1.0]
            ],
			c: [
                [Float64.NaN, Float64.NaN, Float64.NaN],
                [Float64.NaN, Float64.NaN, Float64.NaN],
                [Float64.NaN, Float64.NaN, Float64.NaN],
				[Float64.NaN, Float64.NaN, Float64.NaN]
            ],
			alpha: 2.0,
            ans: [
                [134.0, 136.0, 154.0],
                [220.0, 234.0, 266.0],
                [286.0, 308.0, 350.0],
				[326.0, 358.0, 412.0]
			]
		),
		dsymmTestStruct(
            side: Right,
            ul: Upper,
            m: 3,
            n: 4,
			a: [
				[2.0, 0.0, 0.0, 0.0],
				[3.0, 6.0, 0.0, 0.0],
                [4.0, 7.0, 10.0, 0.0],
				[3.0, 4.0, 5.0, 6.0]
            ],
            b: [
				[2.0, 3.0, 4.0, 9.0],
				[5.0, 6.0, 7.0, -3.0],
                [8.0, 9.0, 10.0, -2.0]
            ],
            c: [
				[Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN],
                [Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN],
                [Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN]
			],
			alpha: 2.0,
            ans: [
                [8.0, 36.0, 80.0, 108.0],
                [20.0, 72.0, 140.0, -36.0],
				[32.0, 108.0, 200.0, -24.0]
            ]
		),
		dsymmTestStruct(
            side: Right,
            ul: Upper,
            m: 4,
            n: 3,
			a: [
				[2.0, 0.0, 0.0],
                [3.0, 6.0, 0.0],
                [4.0, 7.0, 10.0]
			],
			b: [
				[2.0, 3.0, 4.0],
                [5.0, 6.0, 7.0],
                [8.0, 9.0, 10.0],
                [2.0, 1.0, 1.0]
            ],
			c: [
				[Float64.NaN, Float64.NaN, Float64.NaN],
				[Float64.NaN, Float64.NaN, Float64.NaN],
				[Float64.NaN, Float64.NaN, Float64.NaN],
				[Float64.NaN, Float64.NaN, Float64.NaN]
			],
			alpha: 2.0,
            ans: [
                [8.0, 36.0, 80.0, 20.0],
				[72.0, 140.0, 32.0, 108.0],
				[200.0, 8.0, 12.0, 20.0]
			]
		),
		dsymmTestStruct(
            side: Right,
            ul: Lower,
            m: 3,
            n: 4,
			a: [
				[2.0, 0.0, 0.0, 0.0],
				[3.0, 6.0, 0.0, 0.0],
                [4.0, 7.0, 10.0, 0.0],
				[3.0, 4.0, 5.0, 6.0]
            ],
            b: [
				[2.0, 3.0, 4.0, 2.0],
				[5.0, 6.0, 7.0, 1.0],
				[8.0, 9.0, 10.0, 1.0]
			],
			c: [
                [Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN],
                [Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN],
                [Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN]
            ],
			alpha: 2.0,
            ans: [
                [70.0, 120.0, 158.0, 100.0],
                [118.0, 208.0, 274.0, 160.0],
                [172.0, 304.0, 400.0, 232.0]
			]
		),
		dsymmTestStruct(
            side: Right,
            ul: Lower,
            m: 4,
            n: 3,
			a: [
				[2.0, 0.0, 0.0],
				[3.0, 6.0, 0.0],
                [4.0, 7.0, 10.0]
            ],
            b: [
                [2.0, 3.0, 4.0],
				[5.0, 6.0, 7.0],
				[8.0, 9.0, 10.0],
				[2.0, 1.0, 1.0]
			],
			c: [
				[Float64.NaN, Float64.NaN, Float64.NaN],
                [Float64.NaN, Float64.NaN, Float64.NaN],
                [Float64.NaN, Float64.NaN, Float64.NaN],
                [Float64.NaN, Float64.NaN, Float64.NaN]
			],
			alpha: 2.0,
            ans: [
                [58.0, 104.0, 138.0],
                [112.0, 200.0, 264.0],
                [166.0, 296.0, 390.0],
				[22.0, 38.0, 50.0]
			]
		)
    ]
    for (i in 0..tests.size) {
        let test = tests[i]
        let aFlat = flatten(test.a)
        let bFlat = flatten(test.b)
        let cFlat = flatten(test.c)
        let ansFlat = flatten(test.ans)
        nFloat64.dsymm(test.side, test.ul, test.m, test.n, test.alpha, aFlat, test.a[0].size, bFlat, test.n, test.beta, cFlat, test.n)
        if (!fEqualApprox(cFlat, ansFlat, 1e-14)) {
            throw Exception("Case ${i}: Want ${ansFlat}, got ${cFlat}.")
        }
    }
    return true
}