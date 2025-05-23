package cjnum.lapack.cjnum

import cjnum.lapack.*

extend LFloat64Implementation {
    /*
    public func dgeqrf(m: Int64, n: Int64, a: Array<Float64>, lda: Int64, tau: Array<Float64>, work: Array<Float64>,
        lwork: Int64): Unit {
    */
    public func dgeqrf(_: Int64, _: Int64, _: Array<Float64>, _: Int64, _: Array<Float64>, _: Array<Float64>,
        _: Int64): Unit {
        // TODO
        throw Exception("not implemented")
    }
}
