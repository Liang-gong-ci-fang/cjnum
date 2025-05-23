package cjnum.lapack.cjnum

import cjnum.lapack.*
import cjnum.blas.*

extend LFloat64Implementation {
    /*
    public func dpocon(uplo: Uplo, n: Int64, a: Array<Float64>, lda: Int64, anorm: Float64, work: Array<Float64>,
        iwork: Array<Int64>): Float64 {
    */
    public func dpocon(_: Uplo, _: Int64, _: Array<Float64>, _: Int64, _: Float64, _: Array<Float64>,
        _: Array<Int64>): Float64 {
        // TODO
        throw Exception("not implemented")
    }
}
