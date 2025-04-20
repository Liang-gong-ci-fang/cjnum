package cjnum.lapack.cjnum

import cjnum.lapack.*
import cjnum.blas.*

extend LFloat64Implementation {
    public func dtrtri(uplo: Uplo, diag: Diag, n: Int64, a: Array<Float64>, lda: Int64): Bool {
        // TODO
        throw Exception("not implemented")
    }
}