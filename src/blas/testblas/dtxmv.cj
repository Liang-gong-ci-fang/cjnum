package cjnum.blas.testblas

import cjnum.blas.*
import cjnum.floats.*

class vec {
    vec(
        let data!: Array<Float64>,
        let inc!: Int64
    ) {}
}

class dtxmvTestStruct {
    dtxmvTestStruct(
        let n!: Int64,
        let k!: Int64,
        let ul!: Uplo,
        let d!: Diag,
        let ldab!: Int64,
        let tr!: Array<Float64>,
        let tb!: Array<Float64>,
        let tp!: Array<Float64>,
        let ins!: Array<vec>,
        let solNoTrans!: Array<Float64>,
        let solTrans!: Array<Float64>
    ) {}
}

let tests = [
    dtxmvTestStruct(
        n: 3,
        k: 1,
        ul: Upper,
        d: NonUnitDiag,
        tr: [1.0, 2.0, 0.0, 0.0, 3.0, 4.0, 0.0, 0.0, 5.0],
        tb: [1.0, 2.0, 3.0, 4.0, 5.0, 0.0],
        ldab: 2,
        tp: [1.0, 2.0, 0.0, 3.0, 4.0, 5.0],
        ins: [
            vec(
                data: [2.0, 3.0, 4.0],
                inc: 1
            ),
            vec(
                data: [2.0, 1.0, 3.0, 1.0, 4.0],
                inc: 2
            ),
            vec(
                data: [4.0, 1.0, 3.0, 1.0, 2.0],
                inc: -2
            )
        ],
        solNoTrans: [8.0, 25.0, 20.0],
        solTrans:   [2.0, 13.0, 32.0]
    ),
    dtxmvTestStruct(
        n: 3,
        k: 1,
        ul: Upper,
        d: UnitDiag,
        tr: [1.0, 2.0, 0.0, 0.0, 3.0, 4.0, 0.0, 0.0, 5.0],
        tb: [1.0, 2.0, 3.0, 4.0, 5.0, 0.0],
        ldab: 2,
        tp: [1.0, 2.0, 0.0, 3.0, 4.0, 5.0],
        ins: [
            vec(
                data: [2.0, 3.0, 4.0],
                inc: 1
            ),
            vec(
                data: [2.0, 1.0, 3.0, 1.0, 4.0],
                inc: 2
            ),
            vec(
                data: [4.0, 1.0, 3.0, 1.0, 2.0],
                inc: -2
            )
        ],
        solNoTrans: [8.0, 19.0, 4.0],
        solTrans:   [2.0, 7.0, 16.0]
    ),
    dtxmvTestStruct(
        n: 3,
        k: 1,
        ul: Lower,
        d: NonUnitDiag,
        tr: [1.0, 0.0, 0.0, 2.0, 3.0, 0.0, 0.0, 4.0, 5.0],
        tb: [0.0, 1.0, 2.0, 3.0, 4.0, 5.0],
        ldab: 2,
        tp: [1.0, 2.0, 3.0, 0.0, 4.0, 5.0],
        ins: [
            vec(
                data: [2.0, 3.0, 4.0],
                inc: 1
            ),
            vec(
                data: [2.0, 1.0, 3.0, 1.0, 4.0],
                inc: 2
            ),
            vec(
                data: [4.0, 1.0, 3.0, 1.0, 2.0],
                inc: -2
            )
        ],
        solNoTrans: [2.0, 13.0, 32.0],
        solTrans:   [8.0, 25.0, 20.0]
    ),
    dtxmvTestStruct(
        n: 3,
        k: 1,
        ul: Lower,
        d: UnitDiag,
        tr: [1.0, 0.0, 0.0, 2.0, 3.0, 0.0, 0.0, 4.0, 5.0],
        tb: [0.0, 1.0, 2.0, 3.0, 4.0, 5.0],
        ldab: 2,
        tp: [1.0, 2.0, 3.0, 0.0, 4.0, 5.0],
        ins: [
            vec(
                data: [2.0, 3.0, 4.0],
                inc: 1
            ),
            vec(
                data: [2.0, 1.0, 3.0, 1.0, 4.0],
                inc: 2
            ),
            vec(
                data: [4.0, 1.0, 3.0, 1.0, 2.0],
                inc: -2
            )
        ],
        solNoTrans: [2.0, 7.0, 16.0],
        solTrans:   [8.0, 19.0, 4.0]
    )
]

public func dtxmvTest(nFloat64: NFloat64): Bool {
    for (nc in 0..tests.size) {
        let c = tests[nc]
        for (nx in 0..c.ins.size) {
            let x = c.ins[nx]
            var in_ = x.data.clone()
            nFloat64.dtrmv(c.ul, NoTrans, c.d, c.n, c.tr, c.n, in_, x.inc)
            if (!dStridedSliceTolEqual(c.n, in_, x.inc, c.solNoTrans, 1)) {
                throw Exception("Wrong Dtrmv result for: NoTrans  in Case: ${nc}, input: ${nx}")
            }

            in_ = x.data.clone()
            nFloat64.dtrmv(c.ul, Trans, c.d, c.n, c.tr, c.n, in_, x.inc)
            if (!dStridedSliceTolEqual(c.n, in_, x.inc, c.solTrans, 1)) {
                throw Exception("Wrong Dtrmv result for: Trans  in Case: ${nc}, input: ${nx}")
            }

            in_ = x.data.clone()
            nFloat64.dtbmv(c.ul, NoTrans, c.d, c.n, c.k, c.tb, c.ldab, in_, x.inc)
            if (!dStridedSliceTolEqual(c.n, in_, x.inc, c.solNoTrans, 1)) {
                throw Exception("Wrong Dtbmv result for: NoTrans  in Case: ${nc}, input: ${nx}")
            }

            in_ = x.data.clone()
            nFloat64.dtbmv(c.ul, Trans, c.d, c.n, c.k, c.tb, c.ldab, in_, x.inc)
            if (!dStridedSliceTolEqual(c.n, in_, x.inc, c.solTrans, 1)) {
                throw Exception("Wrong Dtbmv result for: Trans  in Case: ${nc}, input: ${nx}")
            }

            in_ = x.data.clone()
            nFloat64.dtpmv(c.ul, NoTrans, c.d, c.n, c.tp, in_, x.inc)
            if (!dStridedSliceTolEqual(c.n, in_, x.inc, c.solNoTrans, 1)) {
                throw Exception("Wrong Dtpmv result for:  NoTrans  in Case: ${nc}, input: ${nx}")
            }

            in_ = x.data.clone()
            nFloat64.dtpmv(c.ul, Trans, c.d, c.n, c.tp, in_, x.inc)
            if (!dStridedSliceTolEqual(c.n, in_, x.inc, c.solTrans, 1)) {
                throw Exception("Wrong Dtpmv result for:  Trans  in Case: ${nc}, input: ${nx}")
            }
        }
    }
    return true
}