package cjnum.lapack.cjnum

import cjnum.blas.blas64.*
import cjnum.blas.*
import std.math.*

extend LFloat64Implementation {
    public func dpotf2(uplo: Uplo, n: Int64, a: Array<Float64>, lda: Int64): Bool {
        if (uplo != Upper && uplo != Lower) {
            throw Exception(badUplo)
        }
        if (n < 0) {
            throw Exception(nLT0)
        }
        if (lda < max(1, n)) {
            throw Exception(badLdA)
        }

        if (n == 0) {
            return true
        }

        if (a.size < (n - 1) * lda + n) {
            throw Exception(shortA)
        }

        let nFloat64 = nFloat64Implementation()

        if (uplo == Upper) {
            for (j in 0..n) {
                var ajj = a[j * lda + j]
                if (j != 0) {
                    ajj -= nFloat64.ddot(j, a[j..], lda, a[j..], lda)
                }
                if (ajj <= 0.0 || ajj.isNaN()) {
                    a[j * lda + j] = ajj
                    return false
                }
                ajj = sqrt(ajj)
                a[j * lda + j] = ajj
                if (j < n - 1) {
                    nFloat64.dgemv(Trans, j, n - j - 1, -1.0, a[j + 1..], lda, a[j..], lda, 1.0, a[j * lda + j + 1..], 1
                    )
                    nFloat64.dscal(n - j - 1, 1.0 / ajj, a[j * lda + j + 1..], 1)
                }
            }
            return true
        }
        for (j in 0..n) {
            var ajj = a[j * lda + j]
            if (j != 0) {
                ajj -= nFloat64.ddot(j, a[j * lda..], 1, a[j * lda..], 1)
            }
            if (ajj <= 0.0 || ajj.isNaN()) {
                a[j * lda + j] = ajj
                return false
            }
            ajj = sqrt(ajj)
            a[j * lda + j] = ajj
            if (j < n - 1) {
                nFloat64.dgemv(NoTrans, n - j - 1, j, -1.0, a[(j + 1) * lda..], lda, a[j * lda..], 1, 1.0,
                    a[(j + 1) * lda + j..], lda)
                nFloat64.dscal(n - j - 1, 1.0 / ajj, a[(j + 1) * lda + j..], lda)
            }
        }
        return true
    }
}
