package cjnum.lapack.cjnum

import cjnum.lapack.*
import cjnum.blas.*

extend LFloat64Implementation {
    // public func dpotri(ul: Uplo, n: Int64, a: Array<Float64>, lda: Int64): Bool {
    public func dpotri(_: Uplo, _: Int64, _: Array<Float64>, _: Int64): Bool {
        // TODO
        throw Exception("not implemented")
    }
}