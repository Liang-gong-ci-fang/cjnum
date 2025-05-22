package cjnum.lapack.cjnum

import std.math.*
import cjnum.lapack.*
import cjnum.blas.*

extend LFloat64Implementation {
    public func dlantb(norm: MatrixNorm, uplo: Uplo, diag: Diag, n: Int64, k: Int64, a: Array<Float64>, lda: Int64, work: Array<Float64>): Float64 {
        if (norm != MaxAbs && norm != MaxRowSum && norm != MaxColumnSum && norm != Frobenius) {
            throw Exception(badNorm)
        } 
        if (uplo != Upper && uplo != Lower) {
            throw Exception(badUplo)
        }
        if (n < 0) {
            throw Exception(nLT0)
        }
        if (k < 0) {
            throw Exception(kdLT0)
        }
        if (lda < k + 1) {
            throw Exception(badLdA)
        }

        if (n == 0) {
            return 0.0
        }

        if (a.size < (n - 1) * lda + k + 1) {
            throw Exception(shortAB)
        }
        if (work.size < n && norm == MaxColumnSum) {
            throw Exception(shortWork)
        }

        var value = 0.0
        if (norm == MaxAbs) {
            if ((uplo == Upper)) {
                var jfirst = 0
                if (((diag == UnitDiag))) {
                    value = 1.0
                    jfirst = 1
                }
                for (i in 0..n) {
                    for (aij in a[i * lda + jfirst..i * lda + min(n - i, k + 1)]) {
                        if (aij.isNaN()) {
                            return aij
                        }
                        let aij_ = abs(aij)
                        if (aij_ > value) {
                            value = aij_
                        }
                    }
                }
            } else {
                var jlast = k + 1
                if ((diag == UnitDiag)) {
                    value = 1.0
                    jlast = k
                }
                for (i in 0..n) {
                    for (aij in a[i * lda + max(0, k - i)..i * lda + jlast]) {
                        if (aij.isNaN()) {
                            return Float64.NaN
                        }
                        let aij_ = abs(aij)
                        if (aij_ > value) {
                            value = aij_
                        }
                    }
                }
            }
        } else if (norm == MaxRowSum) {
            var sum  = 0.0
            if (uplo == Upper) {
                var jfirst = 0
                if ((diag == UnitDiag)) {
                    jfirst = 1
                }
                for (i in 0..n) {
                    sum = 0.0
                    if ((diag == UnitDiag)) {
                        sum = 1.0
                    }
                    for (aij in a[i * lda + jfirst..i * lda + min(n - i, k + 1)]) {
                        sum += abs(aij)
                    }
                    if (sum.isNaN()) {
                        return Float64.NaN
                    }
                    if (sum > value) {
                        value = sum
                    }
                }
            } else {
                var jlast = k + 1
                if ((diag == UnitDiag)) {
                    jlast = k
                }
                for (i in 0..n) {
                    sum = 0.0
                    if ((diag == UnitDiag)) {
                        sum = 1.0
                    }
                    for (aij in a[i * lda + max(0, k - i)..i * lda + jlast]) {
                        sum += abs(aij)
                    }
                    if (sum.isNaN()) {
                        return Float64.NaN
                    }
                    if (sum > value) {
                        value = sum
                    }
                }
            }
        } else if (norm == MaxColumnSum) {
            if ((diag == UnitDiag)) {
                for (i in 0..n) {
                    work[i] = 1.0
                }
            } else {
                for (i in 0..n) {
                    work[i] = 0.0
                }
            }
            if (uplo == Upper) {
                var jfirst = 0
                if ((diag == UnitDiag)) {
                    jfirst = 1
                }
                for (i in 0..n) {
                    let atmp = a[i * lda + jfirst..i * lda + min(n - i, k + 1)]
                    for (j in 0..atmp.size) {
                        let aij = atmp[j]
                        work[i + jfirst + j] += abs(aij)
                    }
                }
            } else {
                var jlast = k + 1
                if ((diag == UnitDiag)) {
                    jlast = k
                }
                for (i in 0..n) {
                    let off = max(0, k-i)
                    let atmp = a[i * lda + off..i * lda + jlast]
                    for (j in 0..atmp.size)  {
                        let aij = atmp[j]
                        work[i + j + off - k] += abs(aij)
                    }
                }
            }
            for (i in 0..n) {
                let wi = work[i]
                if (wi.isNaN()) {
                    return Float64.NaN
                }
                if (wi > value) {
                    value = wi
                }
            }
        } else if (norm == Frobenius) {
            var (scale, sum) = (0.0, 0.0)
            if (uplo == Upper) {
                if ((diag == UnitDiag)) {
                    scale = 1.0
                    sum = Float64(n)
                    if (k > 0) {
                        for (i in 0..n - 1) {
                            let ilen = min(n - i - 1, k)
                            (scale, sum) = dlassq(ilen, a[i * lda + 1..], 1, scale, sum)
                        }
                    }
                } else {
                    scale = 0.0
                    sum = 1.0
                    for (i in 0..n) {
                        let ilen = min(n - i, k + 1)
                        (scale, sum) = dlassq(ilen, a[i * lda..], 1, scale, sum)
                    }
                }
            } else if (uplo == Lower) {
                if ((diag == UnitDiag)) {
                    scale = 1.0
                    sum = Float64(n)
                    if (k > 0) {
                        for (i in 1..n) {
                            let ilen = min(i, k)
                            (scale, sum) = dlassq(ilen, a[i * lda + k - ilen..], 1, scale, sum)
                        }
                    }
                } else {
                    scale = 0.0
                    sum = 1.0
                    for (i in 0..n) {
                        let ilen = min(i, k) + 1
                        (scale, sum) = dlassq(ilen, a[i * lda + k + 1 - ilen..], 1, scale, sum)
                    }
                }
            }
            value = scale * sqrt(sum)
        }

        return value
    }
}