package cjnum.lapack.cjnum

import cjnum.lapack.*
import cjnum.blas.*

extend LFloat64Implementation {
    public func dsyev(jobz: EVJob, uplo: Uplo, n: Int64, a: Array<Float64>, lda: Int64, w: Array<Float64>, work: Array<Float64>, lwork: Int64): Bool {
        // TODO
        throw Exception("not implemented")
    }
}