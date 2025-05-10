package cjnum.lapack.cjnum

import std.math.*

extend LFloat64Implementation {
    public func dlas2(f: Float64, g: Float64, h: Float64): (Float64, Float64) {
        let fa = abs(f)
        let ga = abs(g)
        let ha = abs(h)
        let fhmin = min(fa, ha)
        let fhmax = max(fa, ha)
        if (fhmin == 0.0) {
            if (fhmax == 0.0) {
                return (0.0, ga)
            }
            let v = min(fhmax, ga) / max(fhmax, ga)
            return (0.0, max(fhmax, ga) * sqrt(1.0 + v * v))
        }
        if (ga < fhmax) {
            let as_ = 1.0 + fhmin / fhmax
            let at = (fhmax - fhmin) / fhmax
            let au = (ga / fhmax) * (ga / fhmax)
            let c = 2.0 / (sqrt(as_ * as_ + au) + sqrt(at * at + au))
            return (fhmin * c, fhmax / c)
        }
        let au = fhmax / ga
        if (au == 0.0) {
            return (fhmin * fhmax / ga, ga)
        }
        let as_ = 1.0 + fhmin / fhmax
        let at = (fhmax - fhmin) / fhmax
        let c = 1.0 / (sqrt(1.0 + (as_ * au) * (as_ * au)) + sqrt(1.0 + (at * au) * (at * au)))
        return (2.0 * (fhmin * c) * au, ga / (c + c))
    }
}