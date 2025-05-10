package cjnum.lapack.cjnum

import std.math.*
import cjnum.blas.*
import cjnum.blas.blas64.*
import cjnum.lapack.*

extend LFloat64Implementation {
    public func dlaswp(n: Int64, a: Array<Float64>, lda: Int64, k1: Int64, k2: Int64, ipiv: Array<Int64>, incX: Int64): Unit {
        if (n < 0) {
            throw Exception(nLT0)
        }
        if (k1 < 0) {
            throw Exception(badK1)
        }
        if (k2 < k1) {
            throw Exception(badK2)
        }
        if (lda < max(1, n)) {
            throw Exception(badLdA)
        }
        if (a.size < k2 * lda + n) {
            throw Exception(shortA)
        }
        if (ipiv.size != k2 + 1) {
            throw Exception(badLenIpiv)
        }
        if (incX != 1 && incX != -1) {
            throw Exception(absIncNotOne)
        }

        if (n == 0) {
            return 
        }

        let nFloat64 = nFloat64Implementation()

        if (incX == 1) {
            for (k in k1..=k2) {
                if (k == ipiv[k]) {
                    continue
                }
                nFloat64.dswap(n, a[k * lda..], 1, a[ipiv[k] * lda..], 1)
            }
            return 
        }
        var k = k2
        while (k >= k1) {
            if (k == ipiv[k]) {
                k--
                continue
            }
            nFloat64.dswap(n, a[k * lda..], 1, a[ipiv[k] * lda..], 1)
            k--
        }
    }
}
