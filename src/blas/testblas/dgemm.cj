package cjnum.blas.testblas

import cjnum.blas.*

class dgemmTestStruct {
    dgemmTestStruct(
        let m!: Int64,
        let n!: Int64,
        let k!: Int64,
        let alpha!: Float64,
        let beta!: Float64,
        let a!: Array<Array<Float64>>,
        let b!: Array<Array<Float64>>,
        let c!: Array<Array<Float64>>,
        let ans!: Array<Array<Float64>>
    ) {}
}

let degmmTests = [
    dgemmTestStruct(
        m:     4,
        n:     3,
        k:     2,
        alpha: 2.0,
        beta:  0.5,
        a: [
            [1.0, 2.0],
            [4.0, 5.0],
            [7.0, 8.0],
            [10.0, 11.0]
        ],
        b: [
            [1.0, 5.0, 6.0],
            [5.0, -8.0, 8.0]
        ],
        c: [
            [4.0, 8.0, -9.0],
            [12.0, 16.0, -8.0],
            [1.0, 5.0, 15.0],
            [-3.0, -4.0, 7.0]
        ],
        ans: [
            [24.0, -18.0, 39.5],
            [64.0, -32.0, 124.0],
            [94.5, -55.5, 219.5],
            [128.5, -78.0, 299.5]
        ]
    ),
    dgemmTestStruct(
        m:     4,
        n:     2,
        k:     3,
        alpha: 2.0,
        beta:  0.5,
        a: [
            [1.0, 2.0, 3.0],
            [4.0, 5.0, 6.0],
            [7.0, 8.0, 9.0],
            [10.0, 11.0, 12.0]
        ],
        b: [
            [1.0, 5.0],
            [5.0, -8.0],
            [6.0, 2.0]
        ],
        c: [
            [4.0, 8.0],
            [12.0, 16.0],
            [1.0, 5.0],
            [-3.0, -4.0]
        ],
        ans: [
            [60.0, -6.0],
            [136.0, -8.0],
            [202.5, -19.5],
            [272.5, -30.0]
        ]
    ),
    dgemmTestStruct(
        m:     3,
        n:     2,
        k:     4,
        alpha: 2.0,
        beta:  0.5,
        a: [
            [1.0, 2.0, 3.0, 4.0],
            [4.0, 5.0, 6.0, 7.0],
            [8.0, 9.0, 10.0, 11.0]
        ],
        b: [
            [1.0, 5.0],
            [5.0, -8.0],
            [6.0, 2.0],
            [8.0, 10.0]
        ],
        c: [
            [4.0, 8.0],
            [12.0, 16.0],
            [9.0, -10.0]
        ],
        ans: [
            [124.0, 74.0],
            [248.0, 132.0],
            [406.5, 191.0]
        ]
    ),
    dgemmTestStruct(
        m:     3,
        n:     4,
        k:     2,
        alpha: 2.0,
        beta:  0.5,
        a: [
            [1.0, 2.0],
            [4.0, 5.0],
            [8.0, 9.0]
        ],
        b: [
            [1.0, 5.0, 2.0, 1.0],
            [5.0, -8.0, 2.0, 1.0]
        ],
        c: [
            [4.0, 8.0, 2.0, 2.0],
            [12.0, 16.0, 8.0, 9.0],
            [9.0, -10.0, 10.0, 10.0]
        ],
        ans: [
            [24.0, -18.0, 13.0, 7.0],
            [64.0, -32.0, 40.0, 22.5],
            [110.5, -69.0, 73.0, 39.0]
        ]
    ),
    dgemmTestStruct(
        m:     2,
        n:     4,
        k:     3,
        alpha: 2.0,
        beta:  0.5,
        a: [
            [1.0, 2.0, 3.0],
            [4.0, 5.0, 6.0]
        ],
        b: [
            [1.0, 5.0, 8.0, 8.0],
            [5.0, -8.0, 9.0, 10.0],
            [6.0, 2.0, -3.0, 2.0]
        ],
        c: [
            [4.0, 8.0, 7.0, 8.0],
            [12.0, 16.0, -2.0, 6.0]
        ],
        ans: [
            [60.0, -6.0, 37.5, 72.0],
            [136.0, -8.0, 117.0, 191.0]
        ]
    ),
    dgemmTestStruct(
        m:     2,
        n:     3,
        k:     4,
        alpha: 2.0,
        beta:  0.5,
        a: [
            [1.0, 2.0, 3.0, 4.0],
            [4.0, 5.0, 6.0, 7.0]
        ],
        b: [
            [1.0, 5.0, 8.0],
            [5.0, -8.0, 9.0],
            [6.0, 2.0, -3.0],
            [8.0, 10.0, 2.0]
        ],
        c: [
            [4.0, 8.0, 1.0],
            [12.0, 16.0, 6.0]
        ],
        ans: [
            [124.0, 74.0, 50.5],
            [248.0, 132.0, 149.0]
        ]
    ),
    dgemmTestStruct(
        m:     2,
        n:     3,
        k:     4,
        alpha: 2.0,
        beta:  0.0,
        a: [
            [1.0, 2.0, 3.0, 4.0],
            [4.0, 5.0, 6.0, 7.0]
        ],
        b: [
            [1.0, 5.0, 8.0],
            [5.0, -8.0, 9.0],
            [6.0, 2.0, -3.0],
            [8.0, 10.0, 2.0]
        ],
        c: [
            [4.0, 8.0, 1.0],
            [12.0, 16.0, 6.0]
        ],
        ans: [
            [122.0, 70.0, 50.0],
            [242.0, 124.0, 146.0]
        ]
    )
]

func transpose(a: Array<Array<Float64>>): Array<Array<Float64>> {
    let b = Array<Array<Float64>>(a[0].size, { i => 
        Array<Float64>(a.size, { j => 
            a[j][i]
        })
    })
    return b
}

public func dgemmTest(nFloat64: NFloat64): Bool {
    for (i in 0..degmmTests.size) {
        let test = degmmTests[i]

		dgemmcomp(i, "RowMajorNoTrans", nFloat64, NoTrans, NoTrans,
			test.m, test.n, test.k, test.alpha, test.beta, test.a, test.b, test.c, test.ans)

		dgemmcomp(i, "RowMajorTransA", nFloat64, Trans, NoTrans,
			test.m, test.n, test.k, test.alpha, test.beta, transpose(test.a), test.b, test.c, test.ans)

		dgemmcomp(i, "RowMajorTransB", nFloat64, NoTrans, Trans,
			test.m, test.n, test.k, test.alpha, test.beta, test.a, transpose(test.b), test.c, test.ans)

		dgemmcomp(i, "RowMajorTransBoth", nFloat64, Trans, Trans,
			test.m, test.n, test.k, test.alpha, test.beta, transpose(test.a), transpose(test.b), test.c, test.ans)
    }
    return true
}

func dgemmcomp(i: Int64, name: String, nFloat64: NFloat64, tA: Transpose, tB: Transpose, m: Int64, n: Int64, k: Int64,
                alpha: Float64, beta: Float64, a: Array<Array<Float64>>, b: Array<Array<Float64>>, c: Array<Array<Float64>>, ans: Array<Array<Float64>>): Unit {
	let aFlat = flatten(a)
	let aCopy = flatten(a)
	let bFlat = flatten(b)
	let bCopy = flatten(b)
	let cFlat = flatten(c)
	let ansFlat = flatten(ans)
	let lda = a[0].size
	let ldb = b[0].size
	let ldc = c[0].size

    nFloat64.dgemm(tA, tB, m, n, k, alpha, aFlat, lda, bFlat, ldb, beta, cFlat, ldc)

    if (!dSliceEqual(aFlat, aCopy)) {
        throw Exception("Test ${i} case ${name}: a changed during call to Dgemm")
    }

    if (!dSliceEqual(bFlat, bCopy)) {
        throw Exception("Test ${i} case ${name}: b changed during call to Dgemm")
    }

    if (!dSliceTolEqual(ansFlat, cFlat)) {
        throw Exception("Test ${i} case ${name}: answer mismatch. Expected ${ansFlat}, Found ${cFlat}")
    }
}