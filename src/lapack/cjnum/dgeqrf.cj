package cjnum.lapack.cjnum

import cjnum.blas.*
import cjnum.lapack.*
import std.math.*

extend LFloat64Implementation {
    public func dgeqrf(m: Int64, n: Int64, a: Array<Float64>, lda: Int64, tau: Array<Float64>, work: Array<Float64>,
        lwork: Int64): Unit {
        if (m < 0) {
            throw Exception(mLT0)
        } else if (n < 0) {
            throw Exception(nLT0)
        } else if (lda < max(1, n)) {
            throw Exception(badLdA)
        } else if (lwork < max(1, n) && lwork != -1) {
            throw Exception(badLWork)
        } else if (work.size < max(1, lwork)) {
            throw Exception(shortWork)
        }

        let k = min(m, n)
        if (k == 0) {
            work[0] = 1.0
            return
        }

        var nb = ilaenv(1, "DGEQRF", " ", m, n, -1, -1)
        if (lwork == -1) {
            work[0] = Float64(n * nb)
            return
        }

        if (a.size < (m - 1) * lda + n) {
            throw Exception(shortA)
        }
        if (tau.size != k) {
            throw Exception(badLenTau)
        }

        var nbmin = 2
        var nx: Int64 = 0
        var iws = n
        if (1 < nb && nb < k) {
            nx = max(0, ilaenv(3, "DGEQRF", " ", m, n, -1, -1))
            if (k > nx) {
                iws = n * nb
                if (lwork < iws) {
                    nb = lwork / n
                    nbmin = max(2, ilaenv(2, "DGEQRF", " ", m, n, -1, -1))
                }
            }
        }

        var i = 0
        if (nbmin <= nb && nb < k && nx < k) {
            let ldwork = nb
            while (i < k - nx) {
                let ib = min(k - i, nb)
                dgeqr2(m - i, ib, a[i * lda + i..], lda, tau[i..i + ib], work)
                if (i + ib < n) {
                    dlarft(Forward, ColumnWise, m - i, ib, a[i * lda + i..], lda, tau[i..], work, ldwork)
                    dlarfb(Left, Trans, Forward, ColumnWise, m - i, n - i - ib, ib, a[i * lda + i..], lda, work, ldwork,
                        a[i * lda + i + ib..], lda, work[ib * ldwork..], ldwork)
                }
                i += nb
            }
        }
        if (i < k) {
            dgeqr2(m - i, n - i, a[i * lda + i..], lda, tau[i..], work)
        }
        work[0] = Float64(iws)
    }
}
