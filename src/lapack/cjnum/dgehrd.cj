package cjnum.lapack.cjnum

import cjnum.blas.blas64.*
import cjnum.blas.*
import cjnum.lapack.*
import std.math.*

extend LFloat64Implementation {
    public func dgehrd(n: Int64, ilo: Int64, ihi: Int64, a: Array<Float64>, lda: Int64, tau: Array<Float64>,
        work: Array<Float64>, lwork: Int64): Unit {
        if (n < 0) {
            throw Exception(nLT0)
        } else if (ilo < 0 || max(0, n - 1) < ilo) {
            throw Exception(badIlo)
        } else if (ihi < min(ilo, n - 1) || n <= ihi) {
            throw Exception(badIhi)
        } else if (lda < max(1, n)) {
            throw Exception(badLdA)
        } else if (lwork < max(1, n) && lwork != -1) {
            throw Exception(badLWork)
        } else if (work.size < lwork) {
            throw Exception(shortWork)
        }

        if (n == 0) {
            work[0] = 1.0
            return
        }

        const nbmax = 64
        const ldt = nbmax + 1
        const tsize = ldt * nbmax
        var nb = min(nbmax, ilaenv(1, "DGEHRD", " ", n, ilo, ihi, -1))
        let lwkopt = n * nb + tsize
        if (lwork == -1) {
            work[0] = Float64(lwkopt)
            return
        }

        if (a.size < (n - 1) * lda + n) {
            throw Exception(shortA)
        }
        if (tau.size != n - 1) {
            throw Exception(badLenTau)
        }

        for (i in 0..ilo) {
            tau[i] = 0.0
        }
        for (i in ihi..n - 1) {
            tau[i] = 0.0
        }

        let nh = ihi - ilo + 1
        if (nh <= 1) {
            work[0] = 1.0
            return
        }

        var nbmin = 2
        var nx: Int64 = 0
        if (1 < nb && nb < nh) {
            nx = max(nb, ilaenv(3, "DGEHRD", " ", n, ilo, ihi, -1))
            if (nx < nh) {
                if (lwork < n * nb + tsize) {
                    nbmin = max(2, ilaenv(2, "DGEHRD", " ", n, ilo, ihi, -1))
                    if (lwork >= n * nbmin + tsize) {
                        nb = (lwork - tsize) / n
                    } else {
                        nb = 1
                    }
                }
            }
        }
        let ldwork = nb

        var i: Int64
        if (nb < nbmin || nh <= nb) {
            i = ilo
        } else {
            let nFloat64 = nFloat64Implementation()
            let iwt = n * nb
            i = ilo
            while (i < ihi - nx) {
                let ib = min(nb, ihi - i)

                dlahr2(ihi + 1, i + 1, ib, a[i..], lda, tau[i..], work[iwt..], ldt, work, ldwork)

                let ei = a[(i + ib) * lda + i + ib - 1]
                a[(i + ib) * lda + i + ib - 1] = 1.0
                nFloat64.dgemm(NoTrans, Trans, ihi + 1, ihi - i - ib + 1, ib, -1.0, work, ldwork, a[(i + ib) * lda + i..],
                    lda, 1.0, a[i + ib..], lda)
                a[(i + ib) * lda + i + ib - 1] = ei

                nFloat64.dtrmm(Right, Lower, Trans, UnitDiag, i + 1, ib - 1, 1.0, a[(i + 1) * lda + i..], lda, work,
                    ldwork)
                for (j in 0..=ib - 2) {
                    nFloat64.daxpy(i + 1, -1.0, work[j..], ldwork, a[i + j + 1..], lda)
                }

                dlarfb(Left, Trans, Forward, ColumnWise, ihi - i, n - i - ib, ib, a[(i + 1) * lda + i..], lda,
                    work[iwt..], ldt, a[(i + 1) * lda + i + ib..], lda, work, ldwork)
                i += nb
            }
        }
        dgehd2(n, i, ihi, a, lda, tau, work)
        work[0] = Float64(lwkopt)
    }
}
