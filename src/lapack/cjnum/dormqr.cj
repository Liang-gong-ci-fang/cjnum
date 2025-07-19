package cjnum.lapack.cjnum

import cjnum.lapack.*
import cjnum.blas.*
import std.math.*

extend LFloat64Implementation {
    public func dormqr(side: Side, trans: Transpose, m: Int64, n: Int64, k: Int64, a: Array<Float64>, lda: Int64,
        tau: Array<Float64>, c: Array<Float64>, ldc: Int64, work: Array<Float64>, lwork: Int64): Unit {
        let left = side == Left
        var nq = n
        var nw = m
        if (left) {
            nq = m
            nw = n
        }
        if (!left && side != Right) {
            throw Exception(badSide)
        } else if (trans != NoTrans && trans != Trans) {
            throw Exception(badTrans)
        } else if (m < 0) {
            throw Exception(mLT0)
        } else if (n < 0) {
            throw Exception(nLT0)
        } else if (k < 0) {
            throw Exception(kLT0)
        } else if (left && k > m) {
            throw Exception(kGTM)
        } else if (!left && k > n) {
            throw Exception(kGTN)
        } else if (lda < max(1, k)) {
            throw Exception(badLdA)
        } else if (ldc < max(1, n)) {
            throw Exception(badLdC)
        } else if (lwork < max(1, nw) && lwork != -1) {
            throw Exception(badLWork)
        } else if (work.size < max(1, lwork)) {
            throw Exception(shortWork)
        }

        if (m == 0 || n == 0 || k == 0) {
            work[0] = 1.0
            return
        }

        const nbmax = 64
        const ldt   = nbmax
        const tsize = nbmax * ldt
        let opts = side.toString() + trans.toString()
        var nb = min(nbmax, ilaenv(1, "DORMQR", opts, m, n, k, -1))
        let lworkopt = max(1, nw)*nb + tsize
        if (lwork == -1) {
            work[0] = Float64(lworkopt)
            return
        }

        if (a.size < (nq-1)*lda+k) {
            throw Exception(shortA)
        } else if (tau.size != k) {
            throw Exception(badLenTau)
        } else if (c.size < (m-1)*ldc+n) {
            throw Exception(shortC)
        }

        var nbmin = 2
        if (1 < nb && nb < k) {
            if (lwork < nw*nb+tsize) {
                nb = (lwork - tsize) / nw
                nbmin = max(2, ilaenv(2, "DORMQR", opts, m, n, k, -1))
            }
        }

        if (nb < nbmin || k <= nb) {
            dorm2r(side, trans, m, n, k, a, lda, tau, c, ldc, work)
            work[0] = Float64(lworkopt)
            return
        }

        var ldwork  = nb
        var notrans = trans == NoTrans
        if (left && notrans) {
            for (i in (k - 1) / nb * nb..=0: -nb) {
                let ib = min(nb, k-i)
                dlarft(Forward, ColumnWise, m-i, ib,
                    a[i*lda+i..], lda,
                    tau[i..],
                    work[..tsize], ldt)
                dlarfb(side, trans, Forward, ColumnWise, m-i, n, ib,
                    a[i*lda+i..], lda,
                    work[..tsize], ldt,
                    c[i*ldc..], ldc,
                    work[tsize..], ldwork)
            }
        } else if (left && !notrans) {
            for (i in 0..k: nb) {
                let ib = min(nb, k-i)
                dlarft(Forward, ColumnWise, m-i, ib,
                    a[i*lda+i..], lda,
                    tau[i..],
                    work[..tsize], ldt)
                dlarfb(side, trans, Forward, ColumnWise, m-i, n, ib,
                    a[i*lda+i..], lda,
                    work[..tsize], ldt,
                    c[i*ldc..], ldc,
                    work[tsize..], ldwork)
            }
        } else if (!left && notrans) {
            for (i in 0..k: nb) {
                let ib = min(nb, k-i)
                dlarft(Forward, ColumnWise, n-i, ib,
                    a[i*lda+i..], lda,
                    tau[i..],
                    work[..tsize], ldt)
                dlarfb(side, trans, Forward, ColumnWise, m, n-i, ib,
                    a[i*lda+i..], lda,
                    work[..tsize], ldt,
                    c[i..], ldc,
                    work[tsize..], ldwork)
            }
        } else if (!left && !notrans) {
            for (i in ((k - 1) / nb) * nb..=0: -nb) {
                let ib = min(nb, k-i)
                dlarft(Forward, ColumnWise, n-i, ib,
                    a[i*lda+i..], lda,
                    tau[i..],
                    work[..tsize], ldt)
                dlarfb(side, trans, Forward, ColumnWise, m, n-i, ib,
                    a[i*lda+i..], lda,
                    work[..tsize], ldt,
                    c[i..], ldc,
                    work[tsize..], ldwork)
            }
        }
        work[0] = Float64(lworkopt)
    }
}