package cjnum.lapack.cjnum

import cjnum.lapack.*

extend LFloat64Implementation {
    public func dgeev(jobvl: LeftEVJob, jobvr: RightEVJob, n: Int64, a: Array<Float64>, lda: Int64, wr: Array<Float64>,
        wi: Array<Float64>, vl: Array<Float64>, ldvl: Int64, vr: Array<Float64>, ldvr: Int64, work: Array<Float64>,
        lwork: Int64): Int64 {
        // TODO
        throw Exception("not implemented")
    }
}
