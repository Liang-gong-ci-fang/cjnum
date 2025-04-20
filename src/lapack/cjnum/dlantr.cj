package cjnum.lapack.cjnum

import cjnum.lapack.*
import cjnum.blas.*

extend LFloat64Implementation {
    public func dlantr(norm: MatrixNorm, uplo: Uplo, diag: Diag, m: Int64, n: Int64, a: Array<Float64>, lda: Int64, work: Array<Float64>): Float64 {
        // TODO
        throw Exception("not implemented")
    }
}