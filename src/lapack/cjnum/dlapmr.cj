package cjnum.lapack.cjnum

import cjnum.lapack.*

extend LFloat64Implementation {
    public func dlapmr(forward: Bool, m: Int64, n: Int64, x: Array<Float64>, ldx: Int64, k: Array<Int64>): Unit {
        // TODO
        throw Exception("not implemented")
    }
}