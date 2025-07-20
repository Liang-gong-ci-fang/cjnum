package cjnum.lapack.cjnum

import cjnum.blas.*
import cjnum.lapack.*
import std.math.*

extend LFloat64Implementation {
    public func dgelqf(m: Int64, n: Int64, a: Array<Float64>, lda: Int64, tau: Array<Float64>, work: Array<Float64>,
        lwork: Int64): Unit {
        if (m < 0) {
            throw Exception(mLT0)
        }
        if (n < 0) {
            throw Exception(nLT0)
        }
        if (lda < max(1, n)) {
            throw Exception(badLdA)
        }
        if (lwork < max(1, m) && lwork != -1) {
            throw Exception(badLWork)
        }
        if (work.size < max(1, lwork)) {
            throw Exception(shortWork)
        }

        let k = min(m, n)
        if (k == 0) {
            work[0] = 1.0
            return
        }

        var nb = ilaenv(1, "DGELQF", " ", m, n, -1, -1)
        if (lwork == -1) {
            work[0] = Float64(m * nb)
            return
        }

        if (a.size < (m - 1) * lda + n) {
            throw Exception(shortA)
        }
        if (tau.size < k) {
            throw Exception(shortTau)
        }

        var nbmin = 2
        var nx = 0
        var iws = m
        if (1 < nb && nb < k) {
            nx = max(0, ilaenv(3, "DGELQF", " ", m, n, -1, -1))
            if (nx < k) {
                iws = m * nb
                if (lwork < iws) {
                    nb = lwork / m
                    nbmin = max(2, ilaenv(2, "DGELQF", " ", m, n, -1, -1))
                }
            }
        }
        let ldwork = nb
        var i = 0
        if (nbmin <= nb && nb < k && nx < k) {
            while (i < k - nx) {
                let ib = min(k - i, nb)
                dgelq2(ib, n - i, a[i * lda + i..], lda, tau[i..], work)
                if (i + ib < m) {
                    dlarft(Forward, RowWise, n - i, ib, a[i * lda + i..], lda, tau[i..], work, ldwork)
                    dlarfb(Right, NoTrans, Forward, RowWise, m - i - ib, n - i, ib, a[i * lda + i..], lda, work, ldwork,
                        a[(i + ib) * lda + i..], lda, work[ib * ldwork..], ldwork)
                }
                i += nb
            }
        }
        if (i < k) {
            dgelq2(m - i, n - i, a[i * lda + i..], lda, tau[i..], work)
        }
        work[0] = Float64(iws)
    }
}
