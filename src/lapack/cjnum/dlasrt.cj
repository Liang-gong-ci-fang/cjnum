package cjnum.lapack.cjnum

import std.sort.*
import cjnum.lapack.*

extend LFloat64Implementation {
    public func dlasrt(s: Sort, n: Int64, d: Array<Float64>): Unit {
        if (n < 0) {
            throw Exception(nLT0)
        } else if (d.size < n) {
            throw Exception(shortD)
        }

        if (s == SortIncreasing) {
            stableSort(d[..n])
        } else if (s == SortDecreasing) {
            stableSort(d[..n])
            d[..n].reverse()
        } else {
            throw Exception(badSort)
        }
    }
}
