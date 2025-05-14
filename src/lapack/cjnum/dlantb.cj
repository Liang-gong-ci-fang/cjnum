package cjnum.lapack.cjnum

import cjnum.lapack.*
import cjnum.blas.*

extend LFloat64Implementation {
    public func dlantb(norm: MatrixNorm, uplo: Uplo, diag: Diag, n: Int64, k: Int64, a: Array<Float64>, lda: Int64, work: Array<Float64>): Float64 {
        if (norm != MaxAbs && norm != MaxRowSum && norm != MaxColumnSum && norm != Frobenius) {
            throw Exception(badNorm)
        } 
        if (uplo != Upper && uplo != Lower) {
            throw Exception(badUplo)
        }
        if (n < 0) {
            throw Exception(nLT0)
        }
        if (k < 0) {
            throw Exception(kdLT0)
        }
        if (lda < k + 1) {
            throw Exception(badLdA)
        }
        // TODO
        throw Exception("not implemented")
    }
}