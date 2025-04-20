package cjnum.lapack.cjnum

import cjnum.lapack.*

extend LFloat64Implementation {
    public func dlange(norm: MatrixNorm, m: Int64, n: Int64, a: Array<Float64>, lda: Int64, work: Array<Float64>): Float64 {
        // TODO
        throw Exception("not implemented")
    }
}