package cjnum.lapack.cjnum

import cjnum.lapack.*
import cjnum.blas.*
import std.math.*

extend LFloat64Implementation {
    public func dormlq(side: Side, trans: Transpose, m: Int64, n: Int64, k: Int64, a: Array<Float64>, lda: Int64,
        tau: Array<Float64>, c: Array<Float64>, ldc: Int64, work: Array<Float64>, lwork: Int64): Unit {
        let left = side == Left
        var nw = m
        if (left) {
            nw = n
        }
        if (!left && side != Right) {
            throw Exception(badSide)
        } else if (trans != Trans && trans != NoTrans) {
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
        } else if (left && lda < max(1, m)) {
            throw Exception(badLdA)
        } else if (!left && lda < max(1, n)) {
            throw Exception(badLdA)
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
        const ldt = nbmax
        const tsize = nbmax * ldt

        let opts = side.toString() + trans.toString()
        var nb = min(nbmax, ilaenv(1, "DORMLQ", opts, m, n, k, -1))
        let lworkopt = max(1, nw) * nb + tsize
        if (lwork == -1) {
            work[0] = Float64(lworkopt)
            return
        }

        if (left && a.size < (k - 1) * lda + m) {
            throw Exception(shortA)
        } else if (!left && a.size < (k - 1) * lda + n) {
            throw Exception(shortA)
        } else if (tau.size < k) {
            throw Exception(shortTau)
        } else if (c.size < (m - 1) * ldc + n) {
            throw Exception(shortC)
        }

        var nbmin = 2
        if (1 < nb && nb < k) {
            let iws = nw * nb + tsize
            if (lwork < iws) {
                nb = (lwork - tsize) / nw
                nbmin = max(2, ilaenv(2, "DORMLQ", opts, m, n, k, -1))
            }
        }
        if (nb < nbmin || k <= nb) {
            dorml2(side, trans, m, n, k, a, lda, tau, c, ldc, work)
            work[0] = Float64(lworkopt)
            return
        }

        let t = work[..tsize]
        let wrk = work[tsize..]
        let ldwrk = nb

        let notrans = trans == NoTrans
        var transt = NoTrans
        if (notrans) {
            transt = Trans
        }

        if (left && notrans) {
            for (i in 0..k) {
                let ib = min(nb, k - i)
                dlarft(Forward, RowWise, m - i, ib, a[i * lda + i..], lda, tau[i..], t, ldt)
                dlarfb(side, transt, Forward, RowWise, m - i, n, ib, a[i * lda + i..], lda, t, ldt, c[i * ldc..], ldc,
                    wrk, ldwrk)
            }
        } else if (left && !notrans) {
            for (i in ((k - 1) / nb) * nb..=0 : -1) {
                let ib = min(nb, k - i)
                dlarft(Forward, RowWise, m - i, ib, a[i * lda + i..], lda, tau[i..], t, ldt)
                dlarfb(side, transt, Forward, RowWise, m - i, n, ib, a[i * lda + i..], lda, t, ldt, c[i * ldc..], ldc,
                    wrk, ldwrk)
            }
        } else if (!left && notrans) {
            for (i in ((k - 1) / nb) * nb..=0 : -1) {
                let ib = min(nb, k - i)
                dlarft(Forward, RowWise, n - i, ib, a[i * lda + i..], lda, tau[i..], t, ldt)
                dlarfb(side, transt, Forward, RowWise, m, n - i, ib, a[i * lda + i..], lda, t, ldt, c[i..], ldc, wrk,
                    ldwrk)
            }
        } else if (!left && !notrans) {
            for (i in 0..k) {
                let ib = min(nb, k - i)
                dlarft(Forward, RowWise, n - i, ib, a[i * lda + i..], lda, tau[i..], t, ldt)
                dlarfb(side, transt, Forward, RowWise, m, n - i, ib, a[i * lda + i..], lda, t, ldt, c[i..], ldc, wrk,
                    ldwrk)
            }
        }
        work[0] = Float64(lworkopt)
    }
}
