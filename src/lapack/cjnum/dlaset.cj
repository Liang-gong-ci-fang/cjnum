package cjnum.lapack.cjnum

import cjnum.blas.*

extend LFloat64Implementation {
    public func dlaset(uplo: Uplo, m: Int64, n: Int64, alpha: Float64, beta: Float64, a: Array<Float64>, lda: Int64): Unit {
        if (m < 0) {
            throw Exception(mLT0)
        }
        if (n < 0) {
            throw Exception(nLT0)
        }
        if (lda < max(1, n)) {
            throw Exception(badLdA)
        }

        let minmn = min(m, n)
        if (minmn == 0) {
            return 
        } 

        if (a.size < (m - 1) * lda + n) {
            throw Exception(shortA)
        }

        if (uplo == Upper) {
            for (i in 0..m) {
                for (j in i + 1..n) {
                    a[i * lda + j] = alpha
                }
            }
        } else if (uplo == Lower) {
            for (i in 0..m) {
                for (j in 0..min(i, n)) {
                    a[i * lda + j] = alpha
                }
            }
        } else {
            for (i in 0..m) {
                for (j in 0..n) {
                    a[i * lda + j] = alpha
                }
            }
        }

        for (i in 0..minmn) {
            a[i * lda + i] = beta
        }
    }
}