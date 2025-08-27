package cjnum.lapack.cjnum

import cjnum.blas.*

extend LFloat64Implementation {
    public func dormhr(side: Side, trans: Transpose, m: Int64, n: Int64, ilo: Int64, ihi: Int64, a: Array<Float64>,
        lda: Int64, tau: Array<Float64>, c: Array<Float64>, ldc: Int64, work: Array<Float64>, lwork: Int64): Unit {
        var nq = n
        var nw = m
        if (side == Left) {
            nq = m
            nw = n
        }
        if (side != Left && side != Right) {
            throw Exception(badSide)
        }
        if (trans != NoTrans && trans != Trans) {
            throw Exception(badTrans)
        }
        if (m < 0) {
            throw Exception(mLT0)
        }
        if (n < 0) {
            throw Exception(nLT0)
        }
        if (ilo < 0 || max(1, nq) <= ilo) {
            throw Exception(badIlo)
        }
        if (ihi < min(ilo, nq - 1) || nq <= ihi) {
            throw Exception(badIhi)
        }
        if (lda < max(1, nq)) {
            throw Exception(badLdA)
        }
        if (lwork < max(1, nw) && lwork != -1) {
            throw Exception(badLWork)
        }
        if (work.size < max(1, lwork)) {
            throw Exception(shortWork)
        }

        if (m == 0 || n == 0) {
            work[0] = 1.0
            return
        }

        let nh = ihi - ilo
        var nb = 0
        if (side == Left) {
            var opts = "LN"
            if (trans == Trans) {
                opts = "LT"
            }
            nb = ilaenv(1, "DORMQR", opts, nh, n, nh, -1)
        } else {
            var opts = "RN"
            if (trans == Trans) {
                opts = "RT"
            }
            nb = ilaenv(1, "DORMQR", opts, m, nh, nh, -1)
        }
        let lwkopt = max(1, nw) * nb
        if (lwork == -1) {
            work[0] = Float64(lwkopt)
            return
        }

        if (nh == 0) {
            work[0] = 1.0
            return
        }

        if (a.size < (nq - 1) * lda + nq) {
            throw Exception(shortA)
        }
        if (c.size < (m - 1) * ldc + n) {
            throw Exception(shortC)
        }
        if (tau.size != nq - 1) {
            throw Exception(badLenTau)
        }

        if (side == Left) {
            dormqr(side, trans, nh, n, nh, a[(ilo + 1) * lda + ilo..], lda, tau[ilo..ihi], c[(ilo + 1) * ldc..], ldc,
                work, lwork)
        } else {
            dormqr(side, trans, m, nh, nh, a[(ilo + 1) * lda + ilo..], lda, tau[ilo..ihi], c[ilo + 1..], ldc, work,
                lwork)
        }
        work[0] = Float64(lwkopt)
    }
}
