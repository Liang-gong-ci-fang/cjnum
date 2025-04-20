package cjnum.lapack.cjnum

import cjnum.lapack.*
import cjnum.blas.*

extend LFloat64Implementation {
    public func dpstrf(uplo: Uplo, n: Int64, a: Array<Float64>, lda: Int64, piv: Array<Int64>, tol:  Float64, work: Array<Float64>): (Int64, Bool) {
        // TODO
        throw Exception("not implemented")
    }
}