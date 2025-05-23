package cjnum.lapack.cjnum

import cjnum.lapack.*
import cjnum.blas.*

extend LFloat64Implementation {
    /*
    public func dormlq(side: Side, trans: Transpose, m: Int64, n: Int64, k: Int64, a: Array<Float64>, lda: Int64,
        tau: Array<Float64>, c: Array<Float64>, ldc: Int64, work: Array<Float64>, lwork: Int64): Unit {
    */
    public func dormlq(_: Side, _: Transpose, _: Int64, _: Int64, _: Int64, _: Array<Float64>, _: Int64,
        _: Array<Float64>, _: Array<Float64>, _: Int64, _: Array<Float64>, _: Int64): Unit {
        // TODO
        throw Exception("not implemented")
    }
}
