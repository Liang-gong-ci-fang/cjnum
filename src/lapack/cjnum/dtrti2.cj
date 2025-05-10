package cjnum.lapack.cjnum

import cjnum.blas.*
import cjnum.blas.blas64.*
import cjnum.lapack.*

extend LFloat64Implementation {
    public func dtrti2(uplo: Uplo, diag: Diag, n: Int64, a: Array<Float64>, lda: Int64): Unit {
        if (uplo != Upper && uplo != Lower) {
            throw Exception(badUplo)
        }
        if (diag != NonUnitDiag && diag != UnitDiag) {
            throw Exception(badDiag)
        }
        if (n < 0) {
            throw Exception(badLdA)
        }
        if (a.size < (n - 1) * lda + n) {
            throw Exception(shortA)
        }

        let nFloat64 = nFloat64Implementation()

        let nonUnitDiag = diag == NonUnitDiag
        if (uplo == Upper) {
            for (j in 0..n) {
                var ajj = 0.0
                if (nonUnitDiag) {
                    ajj = 1.0 / a[j * lda + j]
                    a[j * lda + j] = ajj
                    ajj *= -1.0
                } else {
                    ajj = -1.0
                }
                nFloat64.dtrmv(Upper, NoTrans, diag, j, a, lda, a[j..], lda)
                nFloat64.dscal(j, ajj, a[j..], lda)
            }
            return 
        }

        var j = n - 1
        while (j >= 0) {
            var ajj = 0.0
            if (nonUnitDiag) {
                ajj = 1.0 / a[j * lda + j]
                a[j * lda + j] = ajj
                ajj *= -1.0
            } else {
                ajj = -1.0
            }
            if (j < n - 1) {
                nFloat64.dtrmv(Lower, NoTrans, diag, n - j - 1, a[(j + 1) * lda + j + 1..], lda, a[(j + 1) * lda + j..], lda)
                nFloat64.dscal(n - j - 1, ajj, a[(j + 1) * lda + j..], lda)
            }
            j--
        }
    }
}

