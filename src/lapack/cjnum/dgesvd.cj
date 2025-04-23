package cjnum.lapack.cjnum

import cjnum.lapack.*

extend LFloat64Implementation {
    public func dgesvd(jobU: SVDJob, jobVT: SVDJob, m: Int64, n: Int64, a: Array<Float64>, lda: Int64, s: Array<Float64>, u: Array<Float64>, ldu: Int64, vt: Array<Float64>, ldvt: Int64, work: Array<Float64>, lwork: Int64): Bool {
        // TODO
        throw Exception("not implemented")
    }
}