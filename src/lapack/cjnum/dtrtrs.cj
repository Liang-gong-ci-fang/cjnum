package cjnum.lapack.cjnum

import std.math.*
import cjnum.lapack.*
import cjnum.blas.*
import cjnum.blas.blas64.*

extend LFloat64Implementation {
    public func dtrtrs(uplo: Uplo, trans: Transpose, diag: Diag, n: Int64, nrhs: Int64, a: Array<Float64>, lda: Int64,
        b: Array<Float64>, ldb: Int64): Bool {
        if (uplo != Upper && uplo != Lower) {
            throw Exception(badUplo)
        }
        if (trans != NoTrans && trans != Trans && trans != ConjTrans) {
            throw Exception(badTrans)
        }
        if (diag != NonUnitDiag && diag != UnitDiag) {
            throw Exception(badDiag)
        }
        if (n < 0) {
            throw Exception(nLT0)
        }
        if (nrhs < 0) {
            throw Exception(nrhsLT0)
        }
        if (lda < max(1, n)) {
            throw Exception(badLdA)
        }
        if (ldb < max(1, nrhs)) {
            throw Exception(badLdB)
        }

        if (n == 0) {
            return true
        }

        if (a.size < (n - 1) * lda + n) {
            throw Exception(shortA)
        }
        if (b.size < (n - 1) * ldb + nrhs) {
            throw Exception(shortB)
        }

        let nounit = (diag == NonUnitDiag)
        if (nounit) {
            for (i in 0..n) {
                if (a[i * lda + i] == 0.0) {
                    return false
                }
            }
        }
        let nFloat64 = nFloat64Implementation()
        nFloat64.dtrsm(Left, uplo, trans, diag, n, nrhs, 1.0, a, lda, b, ldb)
        return true
    }
}
