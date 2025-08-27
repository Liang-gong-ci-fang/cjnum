package cjnum.lapack.cjnum

import cjnum.blas.*

extend LFloat64Implementation {
    public func dorgtr(uplo: Uplo, n: Int64, a: Array<Float64>, lda: Int64, tau: Array<Float64>, work: Array<Float64>,
        lwork: Int64): Unit {
        if (uplo != Upper && uplo != Lower) {
            throw Exception(badUplo)
        } else if (n < 0) {
            throw Exception(nLT0)
        } else if (lda < max(1, n)) {
            throw Exception(badLdA)
        } else if (lwork < max(1, n - 1) && lwork != -1) {
            throw Exception(badLWork)
        } else if (work.size < max(1, lwork)) {
            throw Exception(shortWork)
        }

        if (n == 0) {
            work[0] = 1.0
            return
        }

        var nb: Int64
        if (uplo == Upper) {
            nb = ilaenv(1, "DORGQL", " ", n - 1, n - 1, n - 1, -1)
        } else {
            nb = ilaenv(1, "DORGQR", " ", n - 1, n - 1, n - 1, -1)
        }
        let lworkopt = max(1, n - 1) * nb
        if (lwork == -1) {
            work[0] = Float64(lworkopt)
            return
        }

        if (a.size < (n - 1) * lda + n) {
            throw Exception(shortA)
        } else if (tau.size < n - 1) {
            throw Exception(shortTau)
        }

        if (uplo == Upper) {
            for (j in 0..n - 1) {
                for (i in 0..=j - 1) {
                    a[i * lda + j] = a[i * lda + j + 1]
                }
                a[(n - 1) * lda + j] = 0.0
            }
            for (i in 0..n - 1) {
                a[i * lda + n - 1] = 0.0
            }
            a[(n - 1) * lda + n - 1] = 1.0

            dorgql(n - 1, n - 1, n - 1, a, lda, tau, work, lwork)
        } else {
            for (j in n - 1..0 : -1) {
                a[j] = 0.0
                for (i in j + 1..n) {
                    a[i * lda + j] = a[i * lda + j - 1]
                }
            }
            a[0] = 1.0
            for (i in 1..n) {
                a[i * lda] = 0.0
            }
            if (n > 1) {
                dorgqr(n - 1, n - 1, n - 1, a[lda + 1..], lda, tau[..n - 1], work, lwork)
            }
        }
        work[0] = Float64(lworkopt)
    }
}
