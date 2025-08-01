package cjnum.lapack.cjnum

import std.math.*
import cjnum.lapack.*
import cjnum.util.nanMax

extend LFloat64Implementation {
    public func dlange(norm: MatrixNorm, m: Int64, n: Int64, a: Array<Float64>, lda: Int64, work: Array<Float64>): Float64 {
        if (norm != MaxRowSum && norm != MaxColumnSum && norm != Frobenius && norm != MaxAbs) { 
            throw Exception(badNorm)
        }
        if (m < 0) {
            throw Exception(mLT0)
        }
        if (n < 0) {
            throw Exception(nLT0)
        }
        if (lda < max(1, n)) {
            throw Exception(badLdA)
        }

        if (m == 0 || n == 0) {
            return 0.0
        }

        if (a.size < (m - 1) * lda + n) {
            throw Exception(badLdA)
        }
        if (norm == MaxColumnSum && work.size < n) {
            throw Exception(shortWork)
        }

        if (norm == MaxAbs) {
            var value = 0.0
            for (i in 0..m) {
                for (j in 0..n) {
                    value = nanMax(value, abs(a[i * lda + j]))
                }
            }
            return value
        } else if (norm == MaxColumnSum) {
            for (i in 0..n) {
                work[i] = 0.0
            }
            for (i in 0..m) {
                for (j in 0..n) {
                    work[j] += abs(a[i * lda + j])
                }
            }
            var value = 0.0
            for (i in 0..n) {
                value = nanMax(value, work[i])
            }
            return value
        } else if (norm == MaxRowSum) {
            var value = 0.0
            for (i in 0..m) {
                var sum = 0.0
                for (j in 0..n) {
                    sum += abs(a[i * lda + j])
                }
                value = nanMax(value, sum)
            }
            return value
        } else {
            var scale = 0.0
            var sum = 1.0
            for (i in 0..m) {
                (scale, sum) = dlassq(n, a[i * lda..], 1, scale, sum)
            }
            return scale * sqrt(sum)
        }
    }
}