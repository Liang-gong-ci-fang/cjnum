package cjnum.lapack.cjnum

import cjnum.blas.blas64.*
import cjnum.blas.*
import std.math.*

extend LFloat64Implementation {
    public func dormr2(side: Side, trans: Transpose, m: Int64, n: Int64, k: Int64, a: Array<Float64>, lda: Int64,
        tau: Array<Float64>, c: Array<Float64>, ldc: Int64, work: Array<Float64>): Unit {
        let left = side == Left
        var nq = n
        var nw = m
        if (left) {
            nq = m
            nw = n
        }
        if (!left && side != Right) {
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
        if (k < 0) {
            throw Exception(kLT0)
        }
        if (left && k > m) {
            throw Exception(kGTM)
        }
        if (!left && k > n) {
            throw Exception(kGTN)
        }
        if (lda < max(1, nq)) {
            throw Exception(badLdA)
        }
        if (ldc < max(1, n)) {
            throw Exception(badLdC)
        }

        if (m == 0 || n == 0 || k == 0) {
            return
        }

        if (a.size < (k - 1) * lda + nq) {
            throw Exception(shortA)
        }
        if (tau.size < k) {
            throw Exception(shortTau)
        }
        if (c.size < (m - 1) * ldc + n) {
            throw Exception(shortC)
        }
        if (work.size < nw) {
            throw Exception(shortWork)
        }

        if (left) {
            if (trans == NoTrans) {
                for (i in k - 1..=0 : -1) {
                    let aii = a[i * lda + (m - k + i)]
                    a[i * lda + (m - k + i)] = 1.0
                    dlarf(side, m - k + i + 1, n, a[i * lda..], 1, tau[i], c, ldc, work)
                    a[i * lda + (m - k + i)] = aii
                }
                return
            }
            for (i in 0..k) {
                let aii = a[i * lda + (m - k + i)]
                a[i * lda + (m - k + i)] = 1.0
                dlarf(side, m - k + i + 1, n, a[i * lda..], 1, tau[i], c, ldc, work)
                a[i * lda + (m - k + i)] = aii
            }
            return
        }
        if (trans == NoTrans) {
            for (i in 0..k) {
                let aii = a[i * lda + (n - k + i)]
                a[i * lda + (n - k + i)] = 1.0
                dlarf(side, m, n - k + i + 1, a[i * lda..], 1, tau[i], c, ldc, work)
                a[i * lda + (n - k + i)] = aii
            }
            return
        }
        for (i in k - 1..=0 : -1) {
            let aii = a[i * lda + (n - k + i)]
            a[i * lda + (n - k + i)] = 1.0
            dlarf(side, m, n - k + i + 1, a[i * lda..], 1, tau[i], c, ldc, work)
            a[i * lda + (n - k + i)] = aii
        }
    }
}
