package cjnum.lapack.cjnum

import cjnum.lapack.*

extend LFloat64Implementation {
    public func dggsvd3(jobU: GSVDJob, jobV: GSVDJob, jobQ :GSVDJob, m: Int64, n: Int64, p: Int64, a: Array<Float64>, lda: Int64, b: Array<Float64>, ldb: Int64, alpha: Array<Float64>, beta: Array<Float64>, u: Array<Float64>, ldu: Int64, v: Array<Float64>, ldv: Int64, q: Array<Float64>, ldq: Int64, work: Array<Float64>, lwork: Int64, iwork: Array<Int64>): (Int64, Int64, Bool) {
        // TODO
        throw Exception("not implemented")
    }
}