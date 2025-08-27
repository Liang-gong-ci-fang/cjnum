package cjnum.lapack.cjnum

import cjnum.blas.*
import cjnum.lapack.*

extend LFloat64Implementation {
    public func dorgql(m: Int64, n: Int64, k: Int64, a: Array<Float64>, lda: Int64, tau: Array<Float64>,
        work: Array<Float64>, lwork: Int64): Unit {
        if (m < 0) {
            throw Exception(mLT0)
        } else if (n < 0) {
            throw Exception(nLT0)
        } else if (n > m) {
            throw Exception(nGTM)
        } else if (k < 0) {
            throw Exception(kLT0)
        } else if (k > n) {
            throw Exception(kGTN)
        } else if (lda < max(1, n)) {
            throw Exception(badLdA)
        } else if (lwork < max(1, n) && lwork != -1) {
            throw Exception(badLWork)
        } else if (work.size < max(1, lwork)) {
            throw Exception(shortWork)
        }

        if (n == 0) {
            work[0] = 1.0
            return
        }

        var nb = ilaenv(1, "DORGQL", " ", m, n, k, -1)
        if (lwork == -1) {
            work[0] = Float64(n * nb)
            return
        }

        if (a.size < (m - 1) * lda + n) {
            throw Exception(shortA)
        } else if (tau.size < k) {
            throw Exception(shortTau)
        }

        var nbmin = 2
        var nx = 0
        var ldwork = 0
        var iws = n
        if (1 < nb && nb < k) {
            nx = max(0, ilaenv(3, "DORGQL", " ", m, n, k, -1))
            if (nx < k) {
                iws = n * nb
                if (lwork < iws) {
                    nb = lwork / n
                    nbmin = max(2, ilaenv(2, "DORGQL", " ", m, n, k, -1))
                }
                ldwork = nb
            }
        }

        var kk = 0
        if (nbmin <= nb && nb < k && nx < k) {
            kk = min(k, ((k - nx + nb - 1) / nb) * nb)
            for (i in m - kk..m) {
                for (j in 0..n - kk) {
                    a[i * lda + j] = 0.0
                }
            }
        }

        dorg2l(m - kk, n - kk, k - kk, a, lda, tau, work)
        if (kk > 0) {
            for (i in k - kk..k : nb) {
                let ib = min(nb, k - i)
                if (n - k + i > 0) {
                    dlarft(Backward, ColumnWise, m - k + i + ib, ib, a[n - k + i..], lda, tau[i..], work, ldwork)

                    dlarfb(Left, NoTrans, Backward, ColumnWise, m - k + i + ib, n - k + i, ib, a[n - k + i..], lda, work,
                        ldwork, a, lda, work[ib * ldwork..], ldwork)
                }

                dorg2l(m - k + i + ib, ib, ib, a[n - k + i..], lda, tau[i..], work)

                for (j in n - k + i..n - k + i + ib) {
                    for (l in m - k + i + ib..m) {
                        a[l * lda + j] = 0.0
                    }
                }
            }
        }
        work[0] = Float64(iws)
    }
}
