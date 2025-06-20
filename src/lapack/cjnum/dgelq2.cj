package cjnum.lapack.cjnum

import cjnum.blas.*
import std.math.*

extend LFloat64Implementation {
    public func dgelq2(m: Int64, n: Int64, a: Array<Float64>, lda: Int64, tau: Array<Float64>, work: Array<Float64>): Unit {
        if (m < 0) {
            throw Exception(mLT0)
        } else if (n < 0) {
            throw Exception(nLT0)
        } else if (lda < max(1, n)) {
            throw Exception(badLdA)
        }

        let k = min(m, n)
        if (k == 0) {
            return
        }

        if (a.size < (m - 1) * lda + n) {
            throw Exception(shortA)
        } else if (tau.size < k) {
            throw Exception(shortTau)
        } else if (work.size < m) {
            throw Exception(shortWork)
        }

        for (i in 0..k) {
            (a[i * lda + i], tau[i]) = dlarfg(n - i, a[i * lda + i], a[i * lda + min(i + 1, n - 1)..], 1)
            if (i < m - 1) {
                let aii = a[i * lda + i]
                a[i * lda + i] = 1.0
                dlarf(Right, m - i - 1, n - i, a[i * lda + i..], 1, tau[i], a[(i + 1) * lda + i..], lda, work)
                a[i * lda + i] = aii
            }
        }
    }
}
