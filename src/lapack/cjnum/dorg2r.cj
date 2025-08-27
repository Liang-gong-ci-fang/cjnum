package cjnum.lapack.cjnum

import cjnum.blas.blas64.*
import cjnum.blas.*

extend LFloat64Implementation {
    public func dorg2r(m: Int64, n: Int64, k: Int64, a: Array<Float64>, lda: Int64, tau: Array<Float64>,
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
        } else if (tau.size != k) {
            throw Exception(badLenTau)
        } else if (work.size < n) {
            throw Exception(shortWork)
        }

        let nFloat64 = nFloat64Implementation()

        for (l in 0..m) {
            for (j in k..n) {
                a[l * lda + j] = 0.0
            }
        }
        for (j in k..n) {
            a[j * lda + j] = 1.0
        }
        for (i in k - 1..=0 : -1) {
            for (i in 0..work.size) {
                work[i] = 0.0
            }
            if (i < n - 1) {
                a[i * lda + i] = 1.0
                dlarf(Left, m - i, n - i - 1, a[i * lda + i..], lda, tau[i], a[i * lda + i + 1..], lda, work)
            }
            if (i < m - 1) {
                nFloat64.dscal(m - i - 1, -tau[i], a[(i + 1) * lda + i..], lda)
            }
            a[i * lda + i] = 1.0 - tau[i]
            for (l in 0..i) {
                a[l * lda + i] = 0.0
            }
        }
    }
}
