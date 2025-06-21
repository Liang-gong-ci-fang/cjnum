package cjnum.lapack.cjnum

import std.math.*
import cjnum.blas.*
import cjnum.lapack.*

extend LFloat64Implementation {
    public func dgeql2(m: Int64, n: Int64, a: Array<Float64>, lda: Int64, tau: Array<Float64>, work: Array<Float64>): Unit {
        if (m < 0) {
            throw Exception(mLT0)
        }
        if (n < 0) {
            throw Exception(nLT0)
        }
        if (lda < max(1, n)) {
            throw Exception(badLdA)
        }

        let k = min(m, n)
        if (k == 0) {
            return
        }

        if (a.size < (m - 1) * lda + n) {
            throw Exception(shortA)
        }
        if (tau.size < k) {
            throw Exception(shortTau)
        }
        if (work.size < n) {
            throw Exception(shortWork)
        }

        var aii = 0.0
        var i = k - 1
        while (i >= 0) {
            (aii, tau[i]) = dlarfg(m - k + i + 1, a[(m - k + i) * lda + n - k + i], a[n - k + i..], lda)

            a[(m - k + i) * lda + n - k + i] = 1.0
            dlarf(Left, m - k + i + 1, n - k + i, a[n - k + i..], lda, tau[i], a, lda, work)
            a[(m - k + i) * lda + n - k + i] = aii
            i--
        }
    }
}
