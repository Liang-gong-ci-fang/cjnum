package cjnum.lapack.cjnum

import cjnum.lapack.*
import cjnum.blas.*

extend LFloat64Implementation {
    /*
    public func dpstrf(uplo: Uplo, n: Int64, a: Array<Float64>, lda: Int64, piv: Array<Int64>, tol: Float64,
        work: Array<Float64>): (Int64, Bool) {
    */
    public func dpstrf(_: Uplo, _: Int64, _: Array<Float64>, _: Int64, _: Array<Int64>, _: Float64,
        _: Array<Float64>): (Int64, Bool) {
        // TODO
        throw Exception("not implemented")
    }
}
