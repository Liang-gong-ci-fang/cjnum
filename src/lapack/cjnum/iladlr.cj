package cjnum.lapack.cjnum

import std.math.*
import cjnum.lapack.*
import cjnum.blas.*

extend LFloat64Implementation {
    public func iladlr(m: Int64, n: Int64, a: Array<Float64>, lda: Int64): Int64 {
        if (m < 0) {
            throw Exception(mLT0)
        }
        if (n < 0) {
            throw Exception(nLT0)
        }
        if (lda < max(1, n)) {
            throw Exception(badLdA)
        }

        if (n == 0 || m == 0) {
            return -1
        }

        if (a.size < (m - 1) * lda + n) {
            throw Exception(shortA)
        }

        if (a[(m - 1) * lda] != 0.0 || a[(m - 1) * lda + n - 1] != 0.0) {
            return m - 1
        }

        var i = m - 1
        while (i >= 0) {
            for (j in 0..n) {
                if (a[i * lda + j] != 0.0) {
                    return i
                }
            }
            i--
        }
        return -1
    }
}