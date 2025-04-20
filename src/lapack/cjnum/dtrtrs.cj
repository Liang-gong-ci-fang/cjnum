package cjnum.lapack.cjnum

import cjnum.lapack.*
import cjnum.blas.*

extend LFloat64Implementation {
    public func dtrtrs(uplo: Uplo, trans: Transpose, diag: Diag, n: Int64, nrhs: Int64, a: Array<Float64>, lda: Int64, b: Array<Float64>, ldb: Int64): Bool {
        // TODO
        throw Exception("not implemented")
    }
}