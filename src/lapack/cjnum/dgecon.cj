package cjnum.lapack.cjnum

import cjnum.lapack.*

extend LFloat64Implementation <: LFloat64{
    public func dgecon(norm: MatrixNorm, n: Int64, a: Array<Float64>, lda: Int64, anorm: Float64, work: Array<Float64>, iwork: Array<Int64>): Float64 {
        // TODO
        throw Exception("not implemented")
    }
}