package cjnum.lapack.cjnum

import cjnum.blas.*

extend LFloat64Implementation {
    // public func dpotrs(ul: Uplo, n: Int64, nrhs: Int64, a: Array<Float64>, lda: Int64, b: Array<Float64>, ldb: Int64): Unit {
    public func dpotrs(_: Uplo, _: Int64, _: Int64, _: Array<Float64>, _: Int64, _: Array<Float64>, _: Int64): Unit {
        // TODO
        throw Exception("not implemented")
    }
}
