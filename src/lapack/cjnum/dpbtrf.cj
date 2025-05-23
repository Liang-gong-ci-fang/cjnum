package cjnum.lapack.cjnum

import cjnum.lapack.*
import cjnum.blas.*

extend LFloat64Implementation {
    // public func dpbtrf(uplo: Uplo, n: Int64, kd: Int64, ab: Array<Float64>, ldab: Int64): Bool {
    public func dpbtrf(_: Uplo, _: Int64, _: Int64, _: Array<Float64>, _: Int64): Bool {
        // TODO
        throw Exception("not implemented")
    }
}