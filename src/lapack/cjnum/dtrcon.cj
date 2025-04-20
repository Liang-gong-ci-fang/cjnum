package cjnum.lapack.cjnum

import cjnum.lapack.*
import cjnum.blas.*

extend LFloat64Implementation {
    public func dtrcon(norm: MatrixNorm, uplo: Uplo, diag: Diag, n: Int64, a: Array<Float64>, lda: Int64, work: Array<Float64>, iwork: Array<Int64>): Float64 {
        // TODO
        throw Exception("not implemented")
    }
}