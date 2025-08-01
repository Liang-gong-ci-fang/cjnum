package cjnum.lapack.cjnum

import cjnum.internal.asm.f64.*
import cjnum.lapack.*
import std.math.*

extend LFloat64Implementation {
    public func dlangb(norm: MatrixNorm, m: Int64, n: Int64, kl: Int64, ku: Int64, ab: Array<Float64>, ldab: Int64): Float64 {
        let ncol = kl + 1 + ku
        if (norm != MaxAbs && norm != MaxRowSum && norm != MaxColumnSum && norm != Frobenius) {
            throw Exception(badNorm)
        } else if (m < 0) {
            throw Exception(mLT0)
        } else if (n < 0) {
            throw Exception(nLT0)
        } else if (kl < 0) {
            throw Exception(klLT0)
        } else if (ku < 0) {
            throw Exception(kuLT0)
        } else if (ldab < ncol) {
            throw Exception(badLdA)
        }

        if (m == 0 || n == 0) {
            return 0.0
        }

        if (ab.size < min(m, n + kl) * ldab) {
            throw Exception(shortAB)
        }

        var value: Float64 = 0.0
        if (norm == MaxAbs) {
            for (i in 0..min(m, n + kl)) {
                let l = max(0, kl - i)
                let u = min(n + kl - i, ncol)
                let abtmp = ab[i * ldab + l..i * ldab + u]
                for (i in 0..abtmp.size) {
                    let aij = abs(abtmp[i])
                    if (aij > value || aij.isNaN()) {
                        value = aij
                    }
                }
            }
        } else if (norm == MaxRowSum) {
            for (i in 0..min(m, n + kl)) {
                let l = max(0, kl - i)
                let u = min(n + kl - i, ncol)
                let sum = l1Norm(ab[i * ldab + l..i * ldab + u])
                if (sum > value || sum.isNaN()) {
                    value = sum
                }
            }
        } else if (norm == MaxColumnSum) {
            for (j in 0..min(m + ku, n)) {
                let jb = min(kl + j, ncol - 1)
                let ib = max(0, j - ku)
                let jlen = min(j + kl, m - 1) - ib + 1
                let sum = l1NormInc(ab[ib * ldab + jb..], jlen, max(1, ldab - 1))
                if (sum > value || sum.isNaN()) {
                    value = sum
                }
            }
        } else if (norm == Frobenius) {
            var scale = 0.0
            var sum = 1.0
            for (i in 0..min(m, n + kl)) {
                let l = max(0, kl - i)
                let u = min(n + kl - i, ncol)
                let ilen = u - l
                (scale, sum) = dlassq(ilen, ab[i * ldab + l..], 1, scale, sum)
            }
            value = scale * sqrt(sum)
        }
        return value
    }
}
