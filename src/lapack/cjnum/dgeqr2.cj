package cjnum.lapack.cjnum

import cjnum.blas.*

extend LFloat64Implementation {
    public func dgeqr2(m: Int64, n: Int64, a: Array<Float64>, lda: Int64, tau: Array<Float64>, work: Array<Float64>): Unit {
        if (m < 0) {
            throw Exception(mLT0)
        } else if (n < 0) {
            throw Exception(nLT0)
        } else if (lda < max(1, n)) {
            throw Exception(badLdA)
        } else if (work.size < n) {
            throw Exception(shortWork)
        }

        let k = min(m, n)
        if (k == 0) {
            return
        }

        if (a.size < (m - 1) * lda + n) {
            throw Exception(shortA)
        } else if (tau.size != k) {
            throw Exception(badLenTau)
        }

        for (i in 0..k) {
            (a[i * lda + i], tau[i]) = dlarfg(m - i, a[i * lda + i], a[min((i + 1), m - 1) * lda + i..], lda)
            if (i < n - 1) {
                let aii = a[i * lda + i]
                a[i * lda + i] = 1.0
                dlarf(Left, m - i, n - i - 1, a[i * lda + i..], lda, tau[i], a[i * lda + i + 1..], lda, work)
                a[i * lda + i] = aii
            }
        }
    }
}
