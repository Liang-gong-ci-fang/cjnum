package cjnum.lapack.testlapack

import std.math.*
import std.random.*
import cjnum.lapack.*

public func dlanhsTest(lFloat64: LFloat64): Bool {
    const tol = 1e-15
    let rnd = Random(1)
    for (n in [0, 1, 2, 4, 9]) {
        for (lda in [max(1, n), n + 5]) {
            let a = randomGeneral(n , n, lda, rnd)
            for (norm in [MaxAbs, MaxRowSum, MaxColumnSum, Frobenius]) {
                let work: Array<Float64> = if (norm == MaxColumnSum) {
                    nanSlice(n)
                } else {
                    []
                }

                let got = lFloat64.dlanhs(norm, a.rows, a.data, lda, work)

                for (i in 2..n) {
                    for (j in 0..max(0, i - 1)) {
                        a.data[i * a.stride + j] = 0.0
                    }
                }
                
                let want = dlange(norm, a.rows, a.cols, a.data, a.stride)

                if (abs(want - got) > tol * want) {
                    throw Exception("Case n=${n},lda=${lda},norm=${normToString(norm)}: unexpected result. Want ${want}, got ${got}.")
                }
            }
        }
    }
    return true
}