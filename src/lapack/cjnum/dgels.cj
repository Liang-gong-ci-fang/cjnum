package cjnum.lapack.cjnum

import cjnum.lapack.*
import cjnum.blas.*

extend LFloat64Implementation {
    /*
    public func dgels(trans: Transpose, m: Int64, n: Int64, nrhs: Int64, a: Array<Float64>, lda: Int64,
        b: Array<Float64>, ldb: Int64, work: Array<Float64>, lwork: Int64): Bool {
    */
    public func dgels(_: Transpose, _: Int64, _: Int64, _: Int64, _: Array<Float64>, _: Int64,
        _: Array<Float64>, _: Int64, _: Array<Float64>, _: Int64): Bool {
        // TODO
        throw Exception("not implemented")
    }
}
