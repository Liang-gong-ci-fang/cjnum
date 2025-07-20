package cjnum.lapack.cjnum

import cjnum.blas.*
import cjnum.lapack.*
import std.math.*

extend LFloat64Implementation {
    public func dormbr(vect: ApplyOrtho, side: Side, trans: Transpose, m: Int64, n: Int64, k: Int64, a: Array<Float64>,
        lda: Int64, tau: Array<Float64>, c: Array<Float64>, ldc: Int64, work: Array<Float64>, lwork: Int64): Unit {
        var nq = n
        var nw = m
        if (side == Left) {
            nq = m
            nw = n
        }
        let applyQ = vect == ApplyQ
        if (!applyQ && vect != ApplyP) {
            throw Exception(badApplyOrtho)
        } else if (side != Left && side != Right) {
            throw Exception(badSide)
        } else if (trans != NoTrans && trans != Trans) {
            throw Exception(badTrans)
        } else if (m < 0) {
            throw Exception(mLT0)
        } else if (n < 0) {
            throw Exception(nLT0)
        } else if (k < 0) {
            throw Exception(kLT0)
        } else if (applyQ && lda < max(1, min(nq, k))) {
            throw Exception(badLdA)
        } else if (!applyQ && lda < max(1, nq)) {
            throw Exception(badLdA)
        } else if (ldc < max(1, n)) {
            throw Exception(badLdC)
        } else if (lwork < max(1, nw) && lwork != -1) {
            throw Exception(badLWork)
        } else if (work.size < max(1, lwork)) {
            throw Exception(shortWork)
        }

        if (m == 0 || n == 0) {
            work[0] = 1.0
            return
        }

        var opts: String
        if (side == Left) {
            opts = "L"
        } else {
            opts = "R"
        }
        if (trans == Trans) {
            opts += "T"
        } else {
            opts += "N"
        }
        var nb: Int
        if (applyQ) {
            if (side == Left) {
                nb = ilaenv(1, "DORMQR", opts, m - 1, n, m - 1, -1)
            } else {
                nb = ilaenv(1, "DORMQR", opts, m, n - 1, n - 1, -1)
            }
        } else {
            if (side == Left) {
                nb = ilaenv(1, "DORMLQ", opts, m - 1, n, m - 1, -1)
            } else {
                nb = ilaenv(1, "DORMLQ", opts, m, n - 1, n - 1, -1)
            }
        }
        let lworkopt = max(1, nw) * nb
        if (lwork == -1) {
            work[0] = Float64(lworkopt)
            return
        }

        let minnqk = min(nq, k)
        if (applyQ && a.size < (nq - 1) * lda + minnqk) {
            throw Exception(shortA)
        } else if (!applyQ && a.size < (minnqk - 1) * lda + nq) {
            throw Exception(shortA)
        } else if (tau.size < minnqk) {
            throw Exception(shortTau)
        } else if (c.size < (m - 1) * ldc + n) {
            throw Exception(shortC)
        }

        if (applyQ) {
            if (nq >= k) {
                dormqr(side, trans, m, n, k, a, lda, tau[..k], c, ldc, work, lwork)
            } else if (nq > 1) {
                var mi = m
                var ni = n - 1
                var i1 = 0
                var i2 = 1
                if (side == Left) {
                    mi = m - 1
                    ni = n
                    i1 = 1
                    i2 = 0
                }
                dormqr(side, trans, mi, ni, nq - 1, a[lda..], lda, tau[..nq - 1], c[i1 * ldc + i2..], ldc, work, lwork)
            }
            work[0] = Float64(lworkopt)
            return
        }

        var transt = Trans
        if (trans == Trans) {
            transt = NoTrans
        }

        if (nq > k) {
            dormlq(side, transt, m, n, k, a, lda, tau, c, ldc, work, lwork)
        } else if (nq > 1) {
            var mi = m
            var ni = n - 1
            var i1 = 0
            var i2 = 1
            if (side == Left) {
                mi = m - 1
                ni = n
                i1 = 1
                i2 = 0
            }
            dormlq(side, transt, mi, ni, nq - 1, a[1..], lda, tau, c[i1 * ldc + i2..], ldc, work, lwork)
        }
        work[0] = Float64(lworkopt)
    }
}
