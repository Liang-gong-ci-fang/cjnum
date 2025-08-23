package cjnum.lapack.cjnum

import std.math.*
import cjnum.lapack.*
import cjnum.blas.*
import cjnum.blas.blas64.*

extend LFloat64Implementation {
    public func dgetrs(trans: Transpose, n: Int64, nrhs: Int64, a: Array<Float64>, lda: Int64, ipiv: Array<Int64>, b: Array<Float64>, ldb: Int64): Unit {
    // public func dgetrs(_: Transpose, _: Int64, _: Int64, _: Array<Float64>, _: Int64, _: Array<Int64>, _: Array<Float64>, _: Int64): Unit {
        // // TODO
        // throw Exception("not implemented")
        if (trans != NoTrans && trans != Trans && trans != ConjTrans) {
            throw Exception(badTrans)
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

        if (n == 0 || nrhs == 0) {
            return 
        }

        if (a.size < (n - 1) * lda + n) {
            throw Exception(shortA)
        }
        if (b.size < (n - 1) * ldb + nrhs) {
            throw Exception(shortB)
        }
        if (ipiv.size != n) {
            throw Exception(badLenIpiv)
        }

        let nFloat64 = nFloat64Implementation()

        if (trans == NoTrans) {
            dlaswp(nrhs, b, ldb, 0, n - 1, ipiv, 1)
            nFloat64.dtrsm(Left, Lower, NoTrans, UnitDiag, n, nrhs, 1.0, a, lda, b, ldb)
            nFloat64.dtrsm(Left, Upper, NoTrans, NonUnitDiag, n, nrhs, 1.0, a, lda, b, ldb)
            return 
        }

        nFloat64.dtrsm(Left, Upper, Trans, NonUnitDiag, n, nrhs, 1.0, a, lda, b, ldb)
        nFloat64.dtrsm(Left, Lower, Trans, UnitDiag, n, nrhs, 1.0, a, lda, b, ldb)
        dlaswp(nrhs, b, ldb, 0, n - 1, ipiv, -1)
    }
}