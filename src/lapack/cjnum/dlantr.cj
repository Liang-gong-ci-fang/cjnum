package cjnum.lapack.cjnum

import std.math.*
import cjnum.lapack.*
import cjnum.blas.*

extend LFloat64Implementation {
    public func dlantr(norm: MatrixNorm, uplo: Uplo, diag: Diag, m: Int64, n: Int64, a: Array<Float64>, lda: Int64, work: Array<Float64>): Float64 {
        if (norm != MaxAbs && norm != MaxRowSum && norm != MaxColumnSum && norm != Frobenius) {
            throw Exception(badNorm)
        } 
        if (uplo != Upper && uplo != Lower) {
            throw Exception(badUplo)
        }
        if (diag != UnitDiag && diag != NonUnitDiag) {
            throw Exception(badDiag)
        }
        if (n < 0) {
            throw Exception(nLT0)
        }
        if (m < 0) {
            throw Exception(mLT0)
        }
        if (lda < max(1, n)) {
            throw Exception(badLdA)
        }

        let minmn = min(m , n)
        if (minmn == 0) {
            return 0.0
        }

        if (a.size < (m - 1) * lda + n) {
            throw Exception(shortA)
        }
        if (norm == MaxColumnSum && work.size < n) {
            throw Exception(shortWork)
        }

        if (norm == MaxAbs) {
            if (diag == UnitDiag) {
                var value = 1.0
                if (uplo == Upper) {
                    for (i in 0..m) {
                        for (j in i + 1..n) {
                            let tmp = abs(a[i * lda + j])
                            if (tmp.isNaN()) {
                                return tmp
                            }
                            if (tmp > value) {
                                value = tmp
                            }
                        }
                    }
                    return value
                }
                for (i in 1..m) {
                    for (j in 0..min(i, n)) {
                        let tmp = abs(a[i * lda + j])
                        if (tmp.isNaN()) {
                            return tmp
                        }
                        if (tmp > value) {
                            value = tmp
                        }
                    }
                }
                return value
            }
            var value = 0.0
            if (uplo == Upper) {
                for (i in 0..m) {
                    for (j in i..n) {
                        let tmp = abs(a[i * lda + j])
                        if (tmp.isNaN()) {
                            return tmp
                        }
                        if (tmp > value) {
                            value = tmp
                        }
                    }
                }
                return value
            }
            for (i in 0..m) {
                for (j in 0..=min(i, n - 1)) {
                    let tmp = abs(a[i * lda + j])
                    if (tmp.isNaN()) {
                        return tmp
                    }
                    if (tmp > value) {
                        value = tmp
                    }
                }
            }
            return value
        } else if (norm == MaxColumnSum) {
            if (diag == UnitDiag) {
                for (i in 0..minmn) {
                    work[i] = 1.0
                }
                for (i in minmn..n) {
                    work[i] = 0.0
                }
                if (uplo == Upper) {
                    for (i in 0..m) {
                        for (j in i + 1..n) {
                            work[j] += abs(a[i * lda + j])
                        }
                    }
                } else {
                    for (i in 1..m) {
                        for (j in 0..min(i, n)) {
                            work[j] += abs(a[i * lda + j])
                        }
                    }
                }
            } else {
                for (i in 0..n) {
                    work[i] = 0.0
                }
                if (uplo == Upper) {
                    for (i in 0..m) {
                        for (j in i..n) {
                            work[j] += abs(a[i * lda + j])
                        }
                    }
                } else {
                    for (i in 0..m) {
                        for (j in 0..=min(i, n - 1)) {
                            work[j] += abs(a[i * lda + j])
                        }
                    }
                }
            }
            var max = 0.0
            for (v in work[..n]) {
                if (v.isNaN()) {
                    return Float64.NaN
                }
                if (v > max) {
                    max = v
                }
            }
            return max
        } else if (norm == MaxRowSum) {
            var maxsum = 0.0
            if (diag == UnitDiag) {
                if (uplo == Upper) {
                    for (i in 0..m) {
                        var sum = 0.0
                        if (i < minmn) {
                            sum = 1.0
                        }
                        for (j in i + 1..n) {
                            sum += abs(a[i * lda + j])
                        }
                        if (sum.isNaN()) {
                            return Float64.NaN
                        }
                        if (sum > maxsum) {
                            maxsum = sum
                        }
                    }
                    return maxsum
                } else {
                    for (i in 0..m) {
                        var sum = 0.0
                        if (i < minmn) {
                            sum = 1.0
                        }
                        for (j in 0..min(i, n)) {
                            sum += abs(a[i * lda + j])
                        }
                        if (sum.isNaN()) {
                            return Float64.NaN
                        }
                        if (sum > maxsum) {
                            maxsum = sum
                        }
                    }
                    return maxsum
                }
            } else {
                if (uplo == Upper) {
                    for (i in 0..m) {
                        var sum = 0.0
                        for (j in i..n) {
                            sum += abs(a[i * lda + j])
                        }
                        if (sum.isNaN()) {
                            return sum
                        }
                        if (sum > maxsum) {
                            maxsum = sum
                        }
                    }
                    return maxsum
                } else {
                    for (i in 0..m) {
                        var sum = 0.0
                        for (j in 0..=min(i, n - 1)) {
                            sum += abs(a[i * lda + j])
                        }
                        if (sum.isNaN()) {
                            return sum
                        }
                        if (sum > maxsum) {
                            maxsum = sum
                        }
                    }
                    return maxsum
                }
            }
        } else {
            var (scale, sum) = (0.0, 0.0)
            if (diag == UnitDiag) {
                scale = 1.0
                sum = Float64(min(m, n))
                if (uplo == Upper) {
                    for (i in 0..min(m, n)) {
                        (scale, sum) = dlassq(n-i-1, a[i * lda + i + 1..], 1, scale, sum)
                    }
                } else {
                    for (i in 1..m) {
                        (scale, sum) = dlassq(min(i, n), a[i * lda..], 1, scale, sum)
                    }
                }
            } else {
                scale = 0.0
                sum = 1.0
                if (uplo == Upper) {
                    for (i in 0..min(m, n)) {
                        (scale, sum) = dlassq(n-i, a[i * lda + i..], 1, scale, sum)
                    }
                } else {
                    for (i in 0..m) {
                        (scale, sum) = dlassq(min(i+1, n), a[i * lda..], 1, scale, sum)
                    }
                }
            }
            return scale * sqrt(sum)
        }
    }
}