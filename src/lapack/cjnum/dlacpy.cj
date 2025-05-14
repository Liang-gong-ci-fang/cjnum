package cjnum.lapack.cjnum

import std.math.*
import cjnum.blas.*

extend LFloat64Implementation {
    public func dlacpy(uplo: Uplo, m: Int64, n: Int64, a: Array<Float64>, lda: Int64, b: Array<Float64>, ldb: Int64): Unit {
        if (uplo != Upper && uplo != Lower && uplo != All) {
            throw Exception(badUplo)
        }
        if (m < 0) {
            throw Exception(mLT0)
        }
        if (n < 0) {
            throw Exception(nLT0)
        }
        if (lda < max(1, n)) {
            throw Exception(badLdA)
        }
        if (ldb < max(1, n)) {
            throw Exception(badLdB)
        }

        if (m == 0 || n == 0) {
            return 
        }

        if (a.size < (m - 1) * lda + n) {
            throw Exception(shortA)
        }
        if (b.size < (m - 1) * ldb + n) {
            throw Exception(shortB)
        }

        if (uplo == Upper) {
            for (i in 0..m) {
                for (j in i..n) {
                    b[i * ldb + j] = a[i * lda + j]
                }
            }
        } else if (uplo == Lower) {
            for (i in 0..m) {
                for (j in 0..min(i + 1, n)) {
                    b[i*ldb+j] = a[i*lda+j]
                }
            }
        } else if (uplo == All) {
            for (i in 0..m) {
                for (j in 0..n) {
                    b[i*ldb+j] = a[i*lda+j]
                }
            }
        }
    }
}