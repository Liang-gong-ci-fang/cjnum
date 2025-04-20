package cjnum.lapack.cjnum

import cjnum.lapack.*

extend LFloat64Implementation {
    public func dorglq(m: Int64, n: Int64, k: Int64, a: Array<Float64>, lda: Int64, tau: Int64, work: Array<Float64>, lwork: Int64): Unit {
        // TODO
        throw Exception("not implemented")
    }
}