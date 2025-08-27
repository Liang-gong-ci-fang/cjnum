package cjnum.lapack.cjnum

import cjnum.blas.*

extend LFloat64Implementation {
    public func dgerq2(m: Int64, n: Int64, a: Array<Float64>, lda: Int64, tau: Array<Float64>, work: Array<Float64>): Unit {
        if (m < 0) {
            throw Exception(mLT0)
        } else if (n < 0) {
            throw Exception(nLT0)
        } else if (lda < max(1, n)) {
            throw Exception(badLdA)
        } else if (work.size < m) {
            throw Exception(shortWork)
        }

        let k = min(m, n)
        if (k == 0) {
            return
        }

        if (a.size < (m - 1) * lda + n) {
            throw Exception(shortA)
        } else if (tau.size < k) {
            throw Exception(shortTau)
        }

        for (i in k - 1..=0 : -1) {
            let mki = m - k + i
            let nki = n - k + i
            var aii: Float64
            (aii, tau[i]) = dlarfg(nki + 1, a[mki * lda + nki], a[mki * lda..], 1)

            a[mki * lda + nki] = 1.0
            dlarf(Right, mki, nki + 1, a[mki * lda..], 1, tau[i], a, lda, work)
            a[mki * lda + nki] = aii
        }
    }
}
