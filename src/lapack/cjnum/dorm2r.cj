package cjnum.lapack.cjnum

import cjnum.blas.*

extend LFloat64Implementation {
    public func dorm2r(side: Side, trans: Transpose, m: Int64, n: Int64, k: Int64, a: Array<Float64>, lda: Int64,
        tau: Array<Float64>, c: Array<Float64>, ldc: Int64, work: Array<Float64>): Unit {
        let left = side == Left
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
        } else if (lda < max(1, k)) {
            throw Exception(badLdA)
        } else if (ldc < max(1, n)) {
            throw Exception(badLdC)
        }

        if (m == 0 || n == 0 || k == 0) {
            return
        }

        if (left && a.size < (m - 1) * lda + k) {
            throw Exception(shortA)
        } else if (!left && a.size < (n - 1) * lda + k) {
            throw Exception(shortA)
        } else if (c.size < (m - 1) * ldc + n) {
            throw Exception(shortC)
        } else if (tau.size != k) {
            throw Exception(badLenTau)
        } else if (left && work.size < n) {
            throw Exception(shortWork)
        } else if (!left && work.size < m) {
            throw Exception(shortWork)
        }

        if (left) {
            if (trans == NoTrans) {
                for (i in k - 1..=0 : -1) {
                    let aii = a[i * lda + i]
                    a[i * lda + i] = 1.0
                    dlarf(side, m - i, n, a[i * lda + i..], lda, tau[i], c[i * ldc..], ldc, work)
                    a[i * lda + i] = aii
                }
                return
            }
            for (i in 0..k) {
                let aii = a[i * lda + i]
                a[i * lda + i] = 1.0
                dlarf(side, m - i, n, a[i * lda + i..], lda, tau[i], c[i * ldc..], ldc, work)
                a[i * lda + i] = aii
            }
            return
        }
        if (trans == NoTrans) {
            for (i in 0..k) {
                let aii = a[i * lda + i]
                a[i * lda + i] = 1.0
                dlarf(side, m, n - i, a[i * lda + i..], lda, tau[i], c[i..], ldc, work)
                a[i * lda + i] = aii
            }
            return
        }
        for (i in k - 1..=0 : -1) {
            let aii = a[i * lda + i]
            a[i * lda + i] = 1.0
            dlarf(side, m, n - i, a[i * lda + i..], lda, tau[i], c[i..], ldc, work)
            a[i * lda + i] = aii
        }
    }
}
