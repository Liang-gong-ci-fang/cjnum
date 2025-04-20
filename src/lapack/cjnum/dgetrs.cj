package cjnum.lapack.cjnum

import cjnum.lapack.*
import cjnum.blas.*

extend LFloat64Implementation {
    public func dgetrs(trans: Transpose, n: Int64, nrhs: Int64, a: Array<Float64>, lda: Int64, ipiv: Array<Int64>, b: Array<Float64>, ldb: Int64): Unit {
        // TODO
        throw Exception("not implemented")
    }
}