package cjnum.lapack.cjnum

import cjnum.blas.*

extend LFloat64Implementation {
    public func dgehd2(n: Int64, ilo: Int64, ihi: Int64, a: Array<Float64>, lda: Int64, tau: Array<Float64>,
        work: Array<Float64>): Unit {
        if (n < 0) {
            throw Exception(nLT0)
        } else if (ilo < 0 || max(0, n - 1) < ilo) {
            throw Exception(badIlo)
        } else if (ihi < min(ilo, n - 1) || n <= ihi) {
            throw Exception(badIhi)
        } else if (lda < max(1, n)) {
            throw Exception(badLdA)
        }

        if (n == 0) {
            return
        }

        if (a.size < (n - 1) * lda + n) {
            throw Exception(shortA)
        } else if (tau.size != n - 1) {
            throw Exception(badLenTau)
        } else if (work.size < n) {
            throw Exception(shortWork)
        }

        for (i in ilo..ihi) {
            var aii = 0.0
            (aii, tau[i]) = dlarfg(ihi - i, a[(i + 1) * lda + i], a[min(i + 2, n - 1) * lda + i..], lda)
            a[(i + 1) * lda + i] = 1.0

            dlarf(Right, ihi + 1, ihi - i, a[(i + 1) * lda + i..], lda, tau[i], a[i + 1..], lda, work)
            dlarf(Left, ihi - i, n - i - 1, a[(i + 1) * lda + i..], lda, tau[i], a[(i + 1) * lda + i + 1..], lda, work)
            a[(i + 1) * lda + i] = aii
        }
    }
}
