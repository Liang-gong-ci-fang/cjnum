package cjnum.lapack.cjnum

import std.math.*
import cjnum.lapack.*
import cjnum.blas.blas64.*
import cjnum.util.*

extend LFloat64Implementation {
    public func dlanhs(norm: MatrixNorm, n: Int64, a: Array<Float64>, lda: Int64, work: Array<Float64>): Float64 {
        if (norm != MaxRowSum && norm != MaxAbs && norm != MaxColumnSum && norm != Frobenius) {
            throw Exception(badNorm)
        }
        if (n < 0) {
            throw Exception(nLT0)
        }
        if (lda < max(1, n)) {
            throw Exception(badLdA)
        }

        if (n == 0) {
            return 0.0
        }

        if (a.size < (n - 1) * lda + n) {
            throw Exception(shortA)
        }
        if (norm == MaxColumnSum && work.size < n) {
            throw Exception(shortWork)
        }

        let nFloat64 = nFloat64Implementation()
        var value = 0.0
        if (norm == MaxAbs) {
            for (i in 0..n) {
                let minj = max(0, i - 1)
                for (v in a[i * lda + minj..i * lda + n]) {
                    value = nanMax(value, abs(v))
                }
            }
        } else if (norm == MaxColumnSum) {
            for (i in 0..n) {
                work[i] = 0.0
            }
            for (i in 0..n) {
                for (j in max(0, i - 1)..n) {
                    work[j] += abs(a[i * lda + j])
                }
            }
            for (v in work[..n]) {
                value = nanMax(value, v)
            }
        } else if (norm == MaxRowSum) {
            for (i in 0..n) {
                let minj = max(0, i - 1)
                let sum = nFloat64.dasum(n- minj, a[i * lda + minj..], 1)
                value = nanMax(value, sum)
            }
        } else if (norm == Frobenius) {
            var scale = 0.0
            var sum = 1.0
            for (i in 0..n) {
                let minj = max(0, i - 1)
                (scale, sum) = dlassq(n - minj, a[i * lda + minj..], 1, scale, sum)
            }
            value = scale * sqrt(sum)
        }

        return value
    }
}