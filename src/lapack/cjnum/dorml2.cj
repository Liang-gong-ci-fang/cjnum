package cjnum.lapack.cjnum

import cjnum.blas.*

extend LFloat64Implementation {
    public func dorml2(side: Side, trans: Transpose, m: Int64, n: Int64, k: Int64, a: Array<Float64>, lda: Int64,
        tau: Array<Float64>, c: Array<Float64>, ldc: Int64, work: Array<Float64>): Unit {
        let left = side == Left
        if (!left && side != Right) {
            throw Exception(badSide)
        }
        if (trans != Trans && trans != NoTrans) {
            throw Exception(badTrans)
        }
        if (m < 0) {
            throw Exception(mLT0)
        }
        if (n < 0) {
            throw Exception(nLT0)
        }
        if (k < 0) {
            throw Exception(kLT0)
        }
        if (left && k > m) {
            throw Exception(kGTM)
        }
        if (!left && k > n) {
            throw Exception(kGTN)
        }
        if (left && lda < max(1, m)) {
            throw Exception(badLdA)
        }
        if (!left && lda < max(1, n)) {
            throw Exception(badLdA)
        }

        if (m == 0 || n == 0 || k == 0) {
            return
        }

        if (left && a.size < (k - 1) * lda + m) {
            throw Exception(shortA)
        }
        if (!left && a.size < (k - 1) * lda + n) {
            throw Exception(shortA)
        }
        if (tau.size < k) {
            throw Exception(shortTau)
        }
        if (c.size < (m - 1) * ldc + n) {
            throw Exception(shortC)
        }
        if (left && work.size < n) {
            throw Exception(shortWork)
        }
        if (!left && work.size < m) {
            throw Exception(shortWork)
        }

        let notrans = trans == NoTrans
        if (left && notrans) {
            for (i in 0..k) {
                let aii = a[i * lda + i]
                a[i * lda + i] = 1.0
                dlarf(side, m - i, n, a[i * lda + i..], 1, tau[i], c[i * ldc..], ldc, work)
                a[i * lda + i] = aii
            }
        } else if (left && !notrans) {
            for (i in k - 1..=0 : -1) {
                let aii = a[i * lda + i]
                a[i * lda + i] = 1.0
                dlarf(side, m - i, n, a[i * lda + i..], 1, tau[i], c[i * ldc..], ldc, work)
                a[i * lda + i] = aii
            }
        } else if (!left && notrans) {
            for (i in k - 1..=0 : -1) {
                let aii = a[i * lda + i]
                a[i * lda + i] = 1.0
                dlarf(side, m, n - i, a[i * lda + i..], 1, tau[i], c[i..], ldc, work)
                a[i * lda + i] = aii
            }
        } else if (!left && !notrans) {
            for (i in 0..k) {
                let aii = a[i * lda + i]
                a[i * lda + i] = 1.0
                dlarf(side, m, n - i, a[i * lda + i..], 1, tau[i], c[i..], ldc, work)
                a[i * lda + i] = aii
            }
        }
    }
}
