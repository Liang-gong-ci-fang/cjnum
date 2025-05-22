package cjnum.lapack.cjnum

import std.math.*
import cjnum.lapack.*
import cjnum.blas.*

extend LFloat64Implementation {
    public func dlansb(norm: MatrixNorm, uplo: Uplo, n: Int64, kd: Int64, ab: Array<Float64>, ldab: Int64, work: Array<Float64>): Float64 {
        if (norm != MaxAbs && norm != MaxRowSum && norm != MaxColumnSum && norm != Frobenius) {
            throw Exception(badNorm)
        }
        if (uplo != Upper && uplo != Lower) {
            throw Exception(badUplo)
        }
        if (n < 0) {
            throw Exception(nLT0)
        }
        if (kd < 0) {
            throw Exception(kdLT0)
        }
        if (ldab < kd + 1) {
            throw Exception(badLdA)
        }

        if (n == 0) {
            return 0.0
        }

        if (ab.size < (n - 1) * ldab + kd + 1) {
            throw Exception(shortAB)
        }
        if (work.size < n && (norm == MaxColumnSum || norm == MaxRowSum)) {
            throw Exception(shortWork)
        }

        var value = 0.0
        if (norm == MaxAbs) {
            if (uplo == Upper) {
                for (i in 0..n) {
                    for (j in 0..min(n - i, kd + 1)) {
                        let aij = abs(ab[i * ldab + j])
                        if (aij > value || aij.isNaN()) {
                            value = aij
                        }
                    }
                }
            } else {
                for (i in 0..n) {
                    for (j in max(0, kd - i)..kd + 1) {
                        let aij = abs(ab[i * ldab + j])
                        if (aij > value || aij.isNaN()) {
                            value = aij
                        }
                    }
                }
            }
        } else if (norm == MaxColumnSum || norm == MaxRowSum) {
            var sum = 0.0
            if (uplo == Upper) {
                for (i in 0..n) {
                    work[i] = 0.0
                }
                for (i in 0..n) {
                    var sum = work[i] + abs(ab[i * ldab])
                    for (j in i + 1..min(i + kd + 1, n)) {
                        let aij = abs(ab[i * ldab + j - i])
                        sum += aij
                        work[j] += aij
                    }
                    if (sum > value || sum.isNaN()) {
                        value = sum
                    }
                }
            } else {
                for (i in 0..n) {
                    sum = 0.0
                    for (j in max(0, i - kd)..i) {
                        let aij = abs(ab[i * ldab + kd + j - i])
                        sum += aij
                        work[j] += aij
                    }
                    work[i] = sum + abs(ab[i * ldab + kd])
                }
                for (i in 0..n) {
                    if (work[i] > value || work[i].isNaN()) {
                        value = work[i]
                    }
                }
            }
        } else if (norm == Frobenius) {
            var scale = 0.0
            var sum = 1.0
            if (uplo == Upper) {
                if (kd > 0) {
                    for (i in 0..n - 1) {
                        let ilen = min(n - i - 1, kd)
                        (scale, sum) = dlassq(ilen, ab[i * ldab + 1..], 1, scale, sum)
                    }
                    sum *= 2.0
                }
                (scale, sum) = dlassq(n, ab, ldab, scale, sum)
            } else {
                if (kd > 0) {
                    for (i in 1..n) {
                        let ilen = min(i, kd)
                        (scale, sum) = dlassq(ilen, ab[i * ldab + kd - ilen..], 1, scale, sum)
                    }
                    sum *= 2.0
                }
                (scale, sum) = dlassq(n, ab[kd..], ldab, scale, sum)
            }
            value = scale * sqrt(sum)
        }
        return value
    }
}