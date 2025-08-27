package cjnum.lapack.testlapack

import cjnum.blas.blas64.General as General64
import cjnum.blas.blas64.*
import cjnum.lapack.*

class dlapmtTestStruct {
    dlapmtTestStruct(
        let forward!: Bool,
        let k!: Array<Int64>,
        let want!: General64
    ) {}
}

public func dlapmtTest(lFloat64: LFloat64): Bool {
    var ti = 0
    for (test in [
            dlapmtTestStruct(
                forward: true,
                k: [0, 1, 2],
                want: General64(
                    rows: 4,
                    cols: 3,
                    stride: 3,
                    data: [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0]
                )
            ),
            dlapmtTestStruct(
                forward: false,
                k: [0, 1, 2],
                want: General64(
                    rows: 4,
                    cols: 3,
                    stride: 3,
                    data: [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0]
                )
            ),
            dlapmtTestStruct(
                forward: true,
                k: [1, 2, 0],
                want: General64(
                    rows: 4,
                    cols: 3,
                    stride: 3,
                    data: [2.0, 3.0, 1.0, 5.0, 6.0, 4.0, 8.0, 9.0, 7.0, 11.0, 12.0, 10.0]
                )
            ),
            dlapmtTestStruct(
                forward: false,
                k: [1, 2, 0],
                want: General64(
                    rows: 4,
                    cols: 3,
                    stride: 3,
                    data: [3.0, 1.0, 2.0, 6.0, 4.0, 5.0, 9.0, 7.0, 8.0, 12.0, 10.0, 11.0]
                )
            )
        ]) {
        let m = test.want.rows
        let n = test.want.cols
        if (test.k.size != n) {
            throw Exception("bad length of k")
        }

        for (extra in [0, 11]) {
            let x = zeros(m, n, n + extra)
            var c = 1
            for (i in 0..m) {
                for (j in 0..n) {
                    x.data[i * x.stride + j] = Float64(c)
                    c++
                }
            }

            let k = Array<Int64>(test.k.size, repeat: 0)
            copy(k, test.k)

            lFloat64.dlapmt(test.forward, m, n, x.data, x.stride, k)

            let prefix = "Case ${ti} (forward=${test.forward},m=${m},n=${n},extra=${extra})"
            if (!generalOutsideAllNaN(x)) {
                throw Exception("${prefix}: out-of-range write to X")
            }

            if (!equalApproxGeneral(x, test.want, 0.0)) {
                throw Exception("${prefix}: unexpected X\n${x}\n${test.want}")
            }
        }
        ti++
    }
    return true
}
