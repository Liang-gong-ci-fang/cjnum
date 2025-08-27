package cjnum.lapack.cjnum

import cjnum.blas.*

extend LFloat64Implementation {
    public func dgebd2(m: Int64, n: Int64, a: Array<Float64>, lda: Int64, d: Array<Float64>, e: Array<Float64>,
        tauQ: Array<Float64>, tauP: Array<Float64>, work: Array<Float64>): Unit {
        if (m < 0) {
            throw Exception(mLT0)
        }
        if (n < 0) {
            throw Exception(nLT0)
        }
        if (lda < max(1, n)) {
            throw Exception(badLdA)
        }

        let minmn = min(m, n)
        if (minmn == 0) {
            return
        }

        if (d.size < minmn) {
            throw Exception(shortD)
        }
        if (e.size < minmn - 1) {
            throw Exception(shortE)
        }
        if (tauQ.size < minmn) {
            throw Exception(shortTauQ)
        }
        if (tauP.size < minmn) {
            throw Exception(shortTauP)
        }
        if (work.size < max(m, n)) {
            throw Exception(shortWork)
        }

        if (m >= n) {
            for (i in 0..n) {
                (a[i * lda + i], tauQ[i]) = dlarfg(m - i, a[i * lda + i], a[min(i + 1, m - 1) * lda + i..], lda)
                d[i] = a[i * lda + i]
                a[i * lda + i] = 1.0

                if (i < n - 1) {
                    dlarf(Left, m - i, n - i - 1, a[i * lda + i..], lda, tauQ[i], a[i * lda + i + 1..], lda, work)
                }
                a[i * lda + i] = d[i]
                if (i < n - 1) {
                    (a[i * lda + i + 1], tauP[i]) = dlarfg(n - i - 1, a[i * lda + i + 1], a[i * lda + min(i + 2, n - 1)..],
                        1)
                    e[i] = a[i * lda + i + 1]
                    a[i * lda + i + 1] = 1.0
                    dlarf(Right, m - i - 1, n - i - 1, a[i * lda + i + 1..], 1, tauP[i], a[(i + 1) * lda + i + 1..], lda,
                        work)
                    a[i * lda + i + 1] = e[i]
                } else {
                    tauP[i] = 0.0
                }
            }
            return
        }

        for (i in 0..m) {
            (a[i * lda + i], tauP[i]) = dlarfg(n - i, a[i * lda + i], a[i * lda + min(i + 1, n - 1)..], 1)
            d[i] = a[i * lda + i]
            a[i * lda + i] = 1.0
            if (i < m - 1) {
                dlarf(Right, m - i - 1, n - i, a[i * lda + i..], 1, tauP[i], a[(i + 1) * lda + i..], lda, work)
            }
            a[i * lda + i] = d[i]
            if (i < m - 1) {
                (a[(i + 1) * lda + i], tauQ[i]) = dlarfg(m - i - 1, a[(i + 1) * lda + i], a[min(i + 2, m - 1) * lda + i..],
                    lda)
                e[i] = a[(i + 1) * lda + i]
                a[(i + 1) * lda + i] = 1.0
                dlarf(Left, m - i - 1, n - i - 1, a[(i + 1) * lda + i..], lda, tauQ[i], a[(i + 1) * lda + i + 1..], lda,
                    work)
                a[(i + 1) * lda + i] = e[i]
            } else {
                tauQ[i] = 0.0
            }
        }
    }
}
