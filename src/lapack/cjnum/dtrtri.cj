package cjnum.lapack.cjnum

import std.math.*
import cjnum.lapack.*
import cjnum.blas.*
import cjnum.blas.blas64.*

extend LFloat64Implementation {
    public func dtrtri(uplo: Uplo, diag: Diag, n: Int64, a: Array<Float64>, lda: Int64): Bool {
        // TODO
        throw Exception("not implemented")
        if (uplo != Upper && uplo != Lower) {
            throw Exception(badUplo)
        }
        if (diag != NonUnitDiag && diag != UnitDiag) {
            throw Exception(badDiag)
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

        if (diag == NonUnitDiag) {
            for (i in 0..n) {
                if (a[i * lda + i] == 0.0) {
                    return false
                }
            }
        }

        let nFloat64 = nFloat64Implementation()

        let nb = ilaenv(1, "DTRTRI", "UD", n, -1, -1, -1)
        if (nb <= 1 || nb > n) {
            dtrti2(uplo, diag, n, a, lda)
            return true
        }
        if (uplo == Upper) {
            var j = 0
            while (j < n) {
                let jb = min(nb, n - j)
                nFloat64.dtrmm(Left, Upper, NoTrans, diag, j, jb, 1.0, a, lda, a[j..], lda)
                nFloat64.dtrsm(Right, Upper, NoTrans, diag, j, jb, -1.0, a[j * lda + j..], lda, a[j..], lda)
                dtrti2(Upper, diag, jb, a[j * lda + j..], lda)
                j += nb
            }
            return true
        }
        let nn = ((n - 1) / nb) * nb
        var j = nn
        while (j >= 0) {
            let jb = min(nb, n - j)
            if (j + jb <= n - 1) {
                nFloat64.dtrmm(Left, Lower, NoTrans, diag, n - j - jb, jb, 1.0, a[(j + jb) * lda + j + jb..], lda, a[(j + jb) * lda + j..], lda)
                nFloat64.dtrsm(Right, Lower, NoTrans, diag, n - j - jb, jb, -1.0, a[j * lda + j..], lda, a[(j + jb) * lda + j..], lda)
            }
            dtrti2(Lower, diag, jb, a[j * lda + j..], lda)
            j -= nb
        }
        return true
    }
}