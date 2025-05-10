package cjnum.lapack.testlapack

import std.random.*
import cjnum.blas.*
import cjnum.lapack.*
import cjnum.blas.blas64.*
import cjnum.blas.blas64.General as General64

class dlaswpTestStruct {
    dlaswpTestStruct(
        let k1!: Int64,
        let k2!: Int64,
        let ipiv!: Array<Int64>,
        let incX!: Int64,
        let want!: General64
    ) {}
}

public func dlaswpTest(lFloat64: LFloat64): Bool {
    let tests = [
        dlaswpTestStruct(
            k1: 0,
            k2: 2,
            ipiv: [0, 1, 2],
            incX: 1,
            want: General64(
                rows: 4,
                cols: 3,
                stride: 3,
                data:  [1.0, 2.0, 3.0,
                        4.0, 5.0, 6.0, 
                        7.0, 8.0, 9.0,
                        10.0, 11.0, 12.0
                ]
            )
        ),
        dlaswpTestStruct(
            k1: 0,
            k2: 2,
            ipiv: [0, 1, 2],
            incX: -1,
            want: General64(
                rows: 4,
                cols: 3,
                stride: 3,
                data:  [1.0, 2.0, 3.0,
                        4.0, 5.0, 6.0, 
                        7.0, 8.0, 9.0,
                        10.0, 11.0, 12.0
                ]
            )
        ),
        dlaswpTestStruct(
            k1: 0,
            k2: 2,
            ipiv: [1, 2, 3],
            incX: 1,
            want: General64(
                rows: 5,
                cols: 3,
                stride: 3,
                data:  [4.0, 5.0, 6.0,
                        7.0, 8.0, 9.0,
                        10.0, 11.0, 12.0,
                        1.0, 2.0, 3.0,
                        13.0, 14.0, 15.0
                ]
            )
        ),
        dlaswpTestStruct(
            k1: 0,
            k2: 2,
            ipiv: [1, 2, 3],
            incX: -1,
            want: General64(
                rows: 5,
                cols: 3,
                stride: 3,
                data:  [10.0, 11.0, 12.0,
                        1.0, 2.0, 3.0,
                        4.0, 5.0, 6.0,
                        7.0, 8.0, 9.0,
                        13.0, 14.0, 15.0
                ]
            )
        )
    ]
    for (i in 0..tests.size) {
        let test = tests[i]
        let m = test.want.rows
        let n = test.want.cols
        let k1 = test.k1
        let k2 = test.k2
        if (test.ipiv.size != k2 + 1) {
            throw Exception("bad length of ipiv")
        }
        let incX = test.incX
        for (extra in [0, 11]) {
            let a = zeros(m, n, n + extra)
            var c = 1
            for (i in 0..m) {
                for (j in 0..n) {
                    a.data[i * a.stride + j] = Float64(c)
                    c++
                }
            }

            let ipiv = test.ipiv.clone()
            lFloat64.dlaswp(n, a.data, a.stride, k1, k2, ipiv, incX)
            
            let prefix = "Case ${i} (m=${m},n=${n},k1=${k1},k2=${k2},extra=${extra})"
            if (!generalOutsideAllNaN(a)) {
                throw Exception("${prefix}: out-of-range write to A")
            }

            if (!equalApproxGeneral(a, test.want, 0.0)) {
                throw Exception("${prefix}: unexpected A")
            }
        }
    }
    return true
}