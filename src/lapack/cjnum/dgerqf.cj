package cjnum.lapack.cjnum

import cjnum.blas.*
import cjnum.lapack.*

extend LFloat64Implementation {
    public func dgerqf(m: Int64, n: Int64, a: Array<Float64>, lda: Int64, tau: Array<Float64>, work: Array<Float64>,
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

        var nb = ilaenv(1, "DGERQF", " ", m, n, -1, -1)
        if (lwork == -1) {
            work[0] = Float64(m * nb)
            return
        }

        if (a.size < (m - 1) * lda + n) {
            throw Exception(shortA)
        }
        if (tau.size != k) {
            throw Exception(badLenTau)
        }

        var nbmin = 2
        var nx = 1
        var iws = m
        var ldwork = 0
        if (1 < nb && nb < k) {
            nx = max(0, ilaenv(3, "DGERQF", " ", m, n, -1, -1))
            if (nx < k) {
                iws = m * nb
                if (lwork < iws) {
                    nb = lwork / m
                    nbmin = max(2, ilaenv(2, "DGERQF", " ", m, n, -1, -1))
                }
                ldwork = nb
            }
        }

        var mu = 0
        var nu = 0
        if (nbmin <= nb && nb < k && nx < k) {
            let ki = ((k - nx - 1) / nb) * nb
            let kk = min(k, ki + nb)

            var i = k - kk + ki
            while (i >= k - kk) {
                let ib = min(k - i, nb)

                dgerq2(ib, n - k + i + ib, a[(m - k + i) * lda..], lda, tau[i..], work)
                if (m - k + i > 0) {
                    dlarft(Backward, RowWise, n - k + i + ib, ib, a[(m - k + i) * lda..], lda, tau[i..], work, ldwork)

                    dlarfb(Right, NoTrans, Backward, RowWise, m - k + i, n - k + i + ib, ib, a[(m - k + i) * lda..], lda,
                        work, ldwork, a, lda, work[ib * ldwork..], ldwork)
                }
                i -= nb
            }
            mu = m - k + i + nb
            nu = n - k + i + nb
        } else {
            mu = m
            nu = n
        }

        if (mu > 0 && nu > 0) {
            dgerq2(mu, nu, a, lda, tau, work)
        }
        work[0] = Float64(iws)
    }
}
