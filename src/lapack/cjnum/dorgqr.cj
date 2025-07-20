package cjnum.lapack.cjnum

import cjnum.blas.*
import cjnum.lapack.*
import std.math.*

extend LFloat64Implementation {
    public func dorgqr(m: Int64, n: Int64, k: Int64, a: Array<Float64>, lda: Int64, tau: Array<Float64>,
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
        } else if (lda < max(1, n) && lwork != -1) {
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

        var nb = ilaenv(1, "DORGQR", " ", m, n, k, -1)

        if (lwork == -1) {
            work[0] = Float64(n * nb)
            return
        }

        if (a.size < (m - 1) * lda + n) {
            throw Exception(shortA)
        } else if (tau.size != k) {
            throw Exception(badLenTau)
        }

        var nbmin = 2
        var nx: Int64 = 0
        let iws = n
        var ldwork: Int64 = 0
        if (1 < nb && nb < k) {
            nx = max(0, ilaenv(3, "DORGQR", " ", m, n, k, -1))
            if (nx < k) {
                ldwork = nb
                let iws = n * ldwork
                if (lwork < iws) {
                    nb = lwork / n
                    ldwork = nb
                    nbmin = max(2, ilaenv(2, "DORGQR", " ", m, n, k, -1))
                }
            }
        }
        var ki: Int64 = 0
        var kk: Int64 = 0
        if (nbmin <= nb && nb < k && nx < k) {
            ki = ((k - nx - 1) / nb) * nb
            kk = min(k, ki + nb)
            for (i in 0..kk) {
                for (j in kk..n) {
                    a[i * lda + j] = 0.0
                }
            }
        }
        if (kk < n) {
            dorg2r(m - kk, n - kk, k - kk, a[kk * lda + kk..], lda, tau[kk..], work)
        }
        if (kk > 0) {
            for (i in ki..=0 : -nb) {
                let ib = min(nb, k - i)
                if (i + ib < n) {
                    dlarft(Forward, ColumnWise, m - i, ib, a[i * lda + i..], lda, tau[i..], work, ldwork)

                    dlarfb(Left, NoTrans, Forward, ColumnWise, m - i, n - i - ib, ib, a[i * lda + i..], lda, work,
                        ldwork, a[i * lda + i + ib..], lda, work[ib * ldwork..], ldwork)
                }
                dorg2r(m - i, ib, ib, a[i * lda + i..], lda, tau[i..i + ib], work)
                for (j in i..i + ib) {
                    for (l in 0..i) {
                        a[l * lda + j] = 0.0
                    }
                }
            }
        }
        work[0] = Float64(iws)
    }
}
