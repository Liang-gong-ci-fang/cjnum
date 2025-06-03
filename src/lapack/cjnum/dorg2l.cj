package cjnum.lapack.cjnum

import cjnum.blas.blas64.*
import cjnum.blas.*
import cjnum.lapack.*
import std.math.*

extend LFloat64Implementation {
    public func dorg2l(m: Int64, n: Int64, k: Int64, a: Array<Float64>, lda: Int64, tau: Array<Float64>,
        work: Array<Float64>): Unit {
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
        } else if (lda < max(1, n)) {
            throw Exception(badLdA)
        }

        if (n == 0) {
            return
        }

        if (a.size < (m - 1) * lda + n) {
            throw Exception(shortA)
        } else if (tau.size < k) {
            throw Exception(shortTau)
        } else if (work.size < n) {
            throw Exception(shortWork)
        }

        for (j in 0..n - k) {
            for (l in 0..m) {
                a[l * lda + j] = 0.0
            }
            a[(m - n + j) * lda + j] = 1.0
        }

        let nFloat64 = nFloat64Implementation()
        for (i in 0..k) {
            let ii = n - k + i

            a[(m - n + ii) * lda + ii] = 1.0
            dlarf(Left, m - n + ii + 1, ii, a[ii..], lda, tau[i], a, lda, work)
            nFloat64.dscal(m - n + ii, -tau[i], a[ii..], lda)
            a[(m - n + ii) * lda + ii] = 1.0 - tau[i]

            for (l in m - n + ii + 1..m) {
                a[l * lda + ii] = 0.0
            }
        }
    }
}
