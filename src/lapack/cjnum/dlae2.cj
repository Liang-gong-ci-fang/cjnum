package cjnum.lapack.cjnum

import std.math.*

extend LFloat64Implementation {
    public func dlae2(a: Float64, b: Float64, c: Float64): (Float64, Float64) {
        let sm = a + c
        let df = a - c
        let adf = abs(df)
        let tb = b + b
        let ab = abs(tb)
        let (acmx, acmn) = if (abs(a) > abs(c)) {
            (a, c)
        } else {
            (c, a)
        }
        let rt = if (adf > ab) {
            adf * sqrt(1.0 + (ab / adf) * (ab / adf))
        } else if (adf < ab) {
            ab * sqrt(1.0 + (adf / ab) * (adf / ab))
        } else {
            ab * sqrt(2.0)
        }
        if (sm < 0.0) {
            let rt1 = 0.5 * (sm - rt)
            let rt2 = (acmx / rt1) * acmn - (b / rt1) * b
            return (rt1, rt2)
        }
        if (sm > 0.0) {
            let rt1 = 0.5 * (sm + rt)
            let rt2 = (acmx / rt1) * acmn - (b / rt1) * b
            return (rt1, rt2)
        }
        let rt1 = 0.5 * rt
        let rt2 = -0.5 * rt
        return (rt1, rt2)
    }
}