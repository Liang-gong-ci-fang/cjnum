package cjnum.lapack.cjnum

import std.math.*

extend LFloat64Implementation {
    public func dlaev2(a: Float64, b: Float64, c: Float64): (Float64, Float64, Float64, Float64) {
        var (rt1, rt2, cs1, sn1) = (0.0, 0.0, 0.0, 0.0)
        let sm = a + c
        let df = a - c
        let adf = abs(df)
        let tb = b + b
        let ab = abs(tb)
        var acmx = c
        var acmn = a
        if (abs(a) > abs(c)) {
            acmx = a
            acmn = c
        }
        var rt = 0.0
        if (adf > ab) {
            rt = adf * sqrt(1.0 + (ab / adf) * (ab / adf))
        } else if (adf < ab) {
            rt = ab * sqrt(1.0 + (adf / ab) * (adf / ab))
        } else {
            rt = ab * sqrt(2.0)
        }
        var sgn1 = 0.0
        if (sm < 0.0) {
            rt1 = 0.5 * (sm - rt)
            sgn1 = -1.0
            rt2 = (acmx / rt1) * acmn - (b / rt1) * b
        } else if (sm > 0.0) {
            rt1 = 0.5 * (sm + rt)
            sgn1 = 1.0
            rt2 = (acmx / rt1) * acmn - (b / rt1) * b
        } else {
            rt1 = 0.5 * rt
            rt2 = -0.5 * rt
            sgn1 = 1.0
        }
        var cs = 0.0
        var sgn2 = 0.0
        if (df >= 0.0) {
            cs = df + rt
            sgn2 = 1.0
        } else {
            cs = df - rt
            sgn2 = -1.0
        }
        let acs = abs(cs)
        if (acs > ab) {
            let ct = -tb / cs
            sn1 = 1.0 / sqrt(1.0 + ct * ct)
            cs1 = ct * sn1
        } else {
            if (ab == 0.0) {
                cs1 = 1.0
                sn1 = 0.0
            } else {
                let tn = -cs / tb
                cs1 = 1.0 / sqrt(1.0 + tn * tn)
                sn1 = tn * cs1
            }
        }
        if (sgn1 == sgn2) {
            (cs1, sn1) = (-sn1, cs1)
        }
        return (rt1, rt2, cs1, sn1)
    }
}
