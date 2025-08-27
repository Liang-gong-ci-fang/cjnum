package cjnum.lapack.cjnum

import cjnum.blas.blas64.*
import cjnum.blas.*

extend LFloat64Implementation {
    public func dlabrd(m: Int64, n: Int64, nb: Int64, a: Array<Float64>, lda: Int64, d: Array<Float64>,
        e: Array<Float64>, tauQ: Array<Float64>, tauP: Array<Float64>, x: Array<Float64>, ldx: Int64, y: Array<Float64>,
        ldy: Int64): Unit {
        if (m < 0) {
            throw Exception(mLT0)
        }
        if (n < 0) {
            throw Exception(nLT0)
        }
        if (nb < 0) {
            throw Exception(nbLT0)
        }
        if (nb > n) {
            throw Exception(nbGTN)
        }
        if (nb > m) {
            throw Exception(nbGTM)
        }
        if (lda < max(1, n)) {
            throw Exception(badLdA)
        }
        if (ldx < max(1, nb)) {
            throw Exception(badLdX)
        }
        if (ldy < max(1, nb)) {
            throw Exception(badLdY)
        }

        if (m == 0 || n == 0 || nb == 0) {
            return
        }

        if (a.size < (m - 1) * lda + n) {
            throw Exception(shortA)
        }
        if (d.size < nb) {
            throw Exception(shortD)
        }
        if (e.size < nb) {
            throw Exception(shortE)
        }
        if (tauQ.size < nb) {
            throw Exception(shortTauQ)
        }
        if (tauP.size < nb) {
            throw Exception(shortTauP)
        }
        if (x.size < (m - 1) * ldx + nb) {
            throw Exception(shortX)
        }
        if (y.size < (n - 1) * ldy + nb) {
            throw Exception(shortY)
        }

        let nFloat64 = nFloat64Implementation()

        if (m >= n) {
            for (i in 0..nb) {
                nFloat64.dgemv(NoTrans, m - i, i, -1.0, a[i * lda..], lda, y[i * ldy..], 1, 1.0, a[i * lda + i..], lda)
                nFloat64.dgemv(NoTrans, m - i, i, -1.0, x[i * ldx..], ldx, a[i..], lda, 1.0, a[i * lda + i..], lda)

                (a[i * lda + i], tauQ[i]) = dlarfg(m - i, a[i * lda + i], a[min(i + 1, m - 1) * lda + i..], lda)
                d[i] = a[i * lda + i]
                if (i < n - 1) {
                    a[i * lda + i] = 1.0
                    nFloat64.dgemv(Trans, m - i, n - i - 1, 1.0, a[i * lda + i + 1..], lda, a[i * lda + i..], lda, 0.0, y[(i + 1) * ldy + i..], ldy)
                    nFloat64.dgemv(Trans, m - i, i, 1.0, a[i * lda..], lda, a[i * lda + i..], lda, 0.0, y[i..], ldy)
                    nFloat64.dgemv(NoTrans, n - i - 1, i, -1.0, y[(i + 1) * ldy..], ldy, y[i..], ldy, 1.0, y[(i + 1) * ldy + i..], ldy)
                    nFloat64.dgemv(Trans, m - i, i, 1.0, x[i * ldx..], ldx, a[i * lda + i..], lda, 0.0, y[i..], ldy)
                    nFloat64.dgemv(Trans, i, n - i - 1, -1.0, a[i + 1..], lda, y[i..], ldy, 1.0, y[(i + 1) * ldy + i..], ldy)
                    nFloat64.dscal(n - i - 1, tauQ[i], y[(i + 1) * ldy + i..], ldy)

                    nFloat64.dgemv(NoTrans, n - i - 1, i + 1, -1.0, y[(i + 1) * ldy..], ldy, a[i * lda..], 1, 1.0, a[i * lda + i + 1..], 1)
                    nFloat64.dgemv(Trans, i, n - i - 1, -1.0, a[i + 1..], lda, x[i * ldx..], 1, 1.0, a[i * lda + i + 1..], 1)

                    (a[i * lda + i + 1], tauP[i]) = dlarfg(n - i - 1, a[i * lda + i + 1], a[i * lda + min(i + 2, n - 1)..], 1)
                    e[i] = a[i * lda + i + 1]
                    a[i * lda + i + 1] = 1.0

                    nFloat64.dgemv(NoTrans, m - i - 1, n - i - 1, 1.0, a[(i + 1) * lda + i + 1..], lda, a[i * lda + i + 1..], 1, 0.0, x[(i + 1) * ldx + i..], ldx)
                    nFloat64.dgemv(Trans, n - i - 1, i + 1, 1.0, y[(i + 1) * ldy..], ldy, a[i * lda + i + 1..], 1, 0.0, x[i..], ldx)
                    nFloat64.dgemv(NoTrans, m - i - 1, i + 1, -1.0, a[(i + 1) * lda..], lda, x[i..], ldx, 1.0, x[(i + 1) * ldx + i..], ldx)
                    nFloat64.dgemv(NoTrans, i, n - i - 1, 1.0, a[i + 1..], lda, a[i * lda + i + 1..], 1, 0.0, x[i..], ldx)
                    nFloat64.dgemv(NoTrans, m - i - 1, i, -1.0, x[(i + 1) * ldx..], ldx, x[i..], ldx, 1.0, x[(i + 1) * ldx + i..], ldx)
                    nFloat64.dscal(m - i - 1, tauP[i], x[(i + 1) * ldx + i..], ldx)
                }
            }
            return
        }

        for (i in 0..nb) {
            nFloat64.dgemv(NoTrans, n - i, i, -1.0, y[i * ldy..], ldy, a[i * lda..], 1, 1.0, a[i * lda + i..], 1)
            nFloat64.dgemv(Trans, i, n - i, -1.0, a[i..], lda, x[i * ldx..], 1, 1.0, a[i * lda + i..], 1)

            (a[i * lda + i], tauP[i]) = dlarfg(n - i, a[i * lda + i], a[i * lda + min(i + 1, n - 1)..], 1)
            d[i] = a[i * lda + i]
            if (i < m - 1) {
                a[i * lda + i] = 1.0
                nFloat64.dgemv(NoTrans, m - i - 1, n - i, 1.0, a[(i + 1) * lda + i..], lda, a[i * lda + i..], 1, 0.0, x[(i + 1) * ldx + i..], ldx)
                nFloat64.dgemv(Trans, n - i, i, 1.0, y[i * ldy..], ldy, a[i * lda + i..], 1, 0.0, x[i..], ldx)
                nFloat64.dgemv(NoTrans, m - i - 1, i, -1.0, a[(i + 1) * lda..], lda, x[i..], ldx, 1.0, x[(i + 1) * ldx + i..], ldx)
                nFloat64.dgemv(NoTrans, i, n - i, 1.0, a[i..], lda, a[i * lda + i..], 1, 0.0, x[i..], ldx)
                nFloat64.dgemv(NoTrans, m - i - 1, i, -1.0, x[(i + 1) * ldx..], ldx, x[i..], ldx, 1.0, x[(i + 1) * ldx + i..], ldx)
                nFloat64.dscal(m - i - 1, tauP[i], x[(i + 1) * ldx + i..], ldx)

                nFloat64.dgemv(NoTrans, m - i - 1, i, -1.0, a[(i + 1) * lda..], lda, y[i * ldy..], 1, 1.0, a[(i + 1) * lda + i..], lda)
                nFloat64.dgemv(NoTrans, m - i - 1, i + 1, -1.0, x[(i + 1) * ldx..], ldx, a[i..], lda, 1.0, a[(i + 1) * lda + i..], lda)

                (a[(i + 1) * lda + i], tauQ[i]) = dlarfg(m - i - 1, a[(i + 1) * lda + i], a[min(i + 2, m - 1) * lda + i..], lda)
                e[i] = a[(i + 1) * lda + i]
                a[(i + 1) * lda + i] = 1.0

                nFloat64.dgemv(Trans, m - i - 1, n - i - 1, 1.0, a[(i + 1) * lda + i + 1..], lda, a[(i + 1) * lda + i..], lda, 0.0, y[(i + 1) * ldy + i..], ldy)
                nFloat64.dgemv(Trans, m - i - 1, i, 1.0, a[(i + 1) * lda..], lda, a[(i + 1) * lda + i..], lda, 0.0, y[i..], ldy)
                nFloat64.dgemv(NoTrans, n - i - 1, i, -1.0, y[(i + 1) * ldy..], ldy, y[i..], ldy, 1.0, y[(i + 1) * ldy + i..], ldy)
                nFloat64.dgemv(Trans, m - i - 1, i + 1, 1.0, x[(i + 1) * ldx..], ldx, a[(i + 1) * lda + i..], lda, 0.0, y[i..], ldy)
                nFloat64.dgemv(Trans, i + 1, n - i - 1, -1.0, a[i + 1..], lda, y[i..], ldy, 1.0, y[(i + 1) * ldy + i..], ldy)
                nFloat64.dscal(n - i - 1, tauQ[i], y[(i + 1) * ldy + i..], ldy)
            }
        }
    }
}
