package cjnum.lapack.cjnum

import cjnum.blas.*

extend LFloat64Implementation {
    /*
    public func dtbtrs(uplo: Uplo, trans: Transpose, diag: Diag, n: Int64, kd: Int64, nrhs: Int64, a: Array<Float64>,
        lda: Int64, b: Array<Float64>, ldb: Int64): Bool {
    */
    public func dtbtrs(_: Uplo, _: Transpose, _: Diag, _: Int64, _: Int64, _: Int64, _: Array<Float64>,
        _: Int64, _: Array<Float64>, _: Int64): Bool {
        // TODO
        throw Exception("not implemented")
    }
}
