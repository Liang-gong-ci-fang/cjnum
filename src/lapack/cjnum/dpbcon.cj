package cjnum.lapack.cjnum

import cjnum.blas.*

extend LFloat64Implementation {
    /*
    public func dpbcon(uplo: Uplo, n: Int64, kd: Int64, ab: Array<Float64>, ldab: Int64, anorm: Float64,
        work: Array<Float64>, iwork: Array<Int64>): Float64 {
    */
    public func dpbcon(_: Uplo, _: Int64, _: Int64, _: Array<Float64>, _: Int64, _: Float64,
        _: Array<Float64>, _: Array<Int64>): Float64 {
        // TODO
        throw Exception("not implemented")
    }
}
