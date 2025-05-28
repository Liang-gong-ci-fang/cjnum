package cjnum.lapack.cjnum

import cjnum.blas.blas64.*
import cjnum.lapack.*
import std.math.*

extend LFloat64Implementation {
    public func dgesc2(n: Int, a: Array<Float64>, lda: Int, rhs: Array<Float64>, ipiv: Array<Int>, jpiv: Array<Int>): Float64 {
        if (n < 0) {
            throw Exception(nLT0)
        } else if (lda < max(1, n)) {
            throw Exception(badLdA)
        }

        if (n == 0) {
            return 0.0
        }

        if (a.size < (n - 1) * lda + n) {
            throw Exception(shortA)
        } else if (rhs.size < n) {
            throw Exception(shortRHS)
        } else if (ipiv.size != n) {
            throw Exception(badLenIpiv)
        } else if (jpiv.size != n) {
            throw Exception(badLenJpiv)
        }

        let smlnum = dlamchS / dlamchP

        dlaswp(1, rhs, 1, 0, n - 1, ipiv[0..n], 1)

        for (i in 0..n - 1) {
            for (j in i + 1..n) {
                rhs[j] -= a[j * lda + i] * rhs[i]
            }
        }

        var scale = 1.0
        let nFloat64 = nFloat64Implementation()
        let i = nFloat64.idamax(n, rhs, 1)
        if (2.0 * smlnum * abs(rhs[i]) > abs(a[(n - 1) * lda + (n - 1)])) {
            let temp = 0.5 / abs(rhs[i])
            nFloat64.dscal(n, temp, rhs, 1)
            scale = scale * temp
        }

        for (i in n - 1..=0 : -1) {
            let temp = 1.0 / a[i * lda + i]
            rhs[i] = rhs[i] * temp
            for (j in i + 1..n) {
                rhs[i] -= rhs[j] * (a[i * lda + j] * temp)
            }
        }

        dlaswp(1, rhs, 1, 0, n - 1, jpiv[0..n], -1)

        return scale
    }
}
