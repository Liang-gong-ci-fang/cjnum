package cjnum.lapack.cjnum

import std.math.*
import cjnum.lapack.*
import cjnum.blas.*

extend LFloat64Implementation {
    public func dlansy(norm: MatrixNorm, uplo: Uplo, n: Int64, a: Array<Float64>, lda: Int64, work: Array<Float64>): Float64 {
        if (norm != MaxRowSum && norm != MaxColumnSum && norm != Frobenius && norm != MaxAbs) {
            throw Exception(badNorm)
        }
        if (uplo != Upper && uplo != Lower) {
            throw Exception(badUplo)
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
        if ((norm == MaxColumnSum || norm == MaxRowSum) && work.size < n) {
            throw Exception(shortWork)
        }

        if (norm == MaxAbs) {
            if (uplo == Upper) {
                var mx = 0.0
                for (i in 0..n) {
                    for (j in i..n) {
                        let v = abs(a[i * lda + j])
                        if (v.isNaN()) {
                            return Float64.NaN
                        }
                        if (v > mx) {
                            mx = v
                        }
                    }
                }
                return mx
            }
            
            var mx = 0.0
            for (i in 0..n) {
                for (j in 0..=i) {
                    let v = abs(a[i * lda + j])
                    if (v.isNaN()) {
                        return Float64.NaN
                    }
                    if (v > mx) {
                        mx = v
                    }
                }
            }
            return mx
        } else if (norm == MaxRowSum || norm == MaxColumnSum) {
            for (i in 0..n) {
                work[i] = 0.0
            }
            if (uplo == Upper) {
                for (i in 0..n) {
                    work[i] += abs(a[i * lda + i])
                    for (j in i + 1..n) {
                        let v = abs(a[i * lda + j])
                        work[i] += v
                        work[j] += v
                    }
                }
            } else {
                for (i in 0..n) {
                    for (j in 0..i) {
                        let v = abs(a[i * lda + j])
                        work[i] += v
                        work[j] += v
                    }
                    work[i] += abs(a[i * lda + i])
                }
            }

            var mx = 0.0
            for (i in 0..n) {
                let v = work[i]
                if (v.isNaN()) {
                    return Float64.NaN
                }
                if (v > mx) {
                    mx = v
                }
            }
            return mx
        } else {
            var scale = 0.0
            var sum = 1.0
            if (uplo == Upper) {
                for (i in 0..n - 1) {
                    (scale, sum) = dlassq(n - i - 1, a[i * lda + i + 1..], 1, scale, sum)
                }
            } else {
                for (i in 1..n) {
                    (scale, sum) = dlassq(i, a[i * lda..], 1, scale, sum)
                }
            }
            sum *= 2.0
            (scale, sum) = dlassq(n, a, lda+1, scale, sum)
            return scale * sqrt(sum)
        }
    }
}