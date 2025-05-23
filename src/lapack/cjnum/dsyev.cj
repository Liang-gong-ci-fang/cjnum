package cjnum.lapack.cjnum

import cjnum.lapack.*
import cjnum.blas.*

extend LFloat64Implementation {
    /*
    public func dsyev(jobz: EVJob, uplo: Uplo, n: Int64, a: Array<Float64>, lda: Int64, w: Array<Float64>,
        work: Array<Float64>, lwork: Int64): Bool {
    */
    public func dsyev(_: EVJob, _: Uplo, _: Int64, _: Array<Float64>, _: Int64, _: Array<Float64>,
        _: Array<Float64>, _: Int64): Bool {
        // TODO
        throw Exception("not implemented")
    }
}
