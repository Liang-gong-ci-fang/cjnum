package cjnum.lapack.cjnum

import cjnum.lapack.*

extend LFloat64Implementation {
    public func dgetri(n: Int64, a: Array<Float64>, lda: Int64, ipiv: Array<Int64>, work: Array<Float64>, lwork: Int64): Bool {
        // TODO
        throw Exception("not implemented")
    }
}