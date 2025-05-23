package cjnum.lapack.cjnum

import cjnum.lapack.*
import cjnum.blas.*

extend LFloat64Implementation {
    /*
    public func dpbtrs(uplo: Uplo, n: Int64, kd: Int64, nrhs: Int64, ab: Array<Float64>, ldab: Int64, b: Array<Float64>,
        ldb: Int64): Unit {
    */
    public func dpbtrs(_: Uplo, _: Int64, _: Int64, _: Int64, _: Array<Float64>, _: Int64, _: Array<Float64>,
        _: Int64): Unit {
        // TODO
        throw Exception("not implemented")
    }
}
