package cjnum.lapack.cjnum

import std.math.*
import cjnum.util.*

extend LFloat64Implementation {
    public func dlartg(f: Float64, g: Float64): (Float64, Float64, Float64) {
        if (g == 0.0) {
            return (1.0, 0.0, f)
        }

        let g1 = abs(g)

        if (f == 0.0) {
            return (0.0, copysign(1.0, g), g1)
        }

        let safmin = dlamchS
        let safmax = 1.0 / safmin
        let rtmin = sqrt(safmin)
        let rtmax = sqrt(safmax / 2.0)

        let f1 = abs(f)

        if (rtmin < f1 && f1 < rtmax && rtmin < g1 && g1 < rtmax) {
            let d = sqrt(f * f + g * g)
            let cs = f1 / d
            let r = copysign(d, f)
            let sn = g / r
            return (cs, sn, r)
        }

        let u = nanMin(nanMax(safmin, nanMax(f1, g1)), safmax)
        let fs = f / u
        let gs = g / u
        let d = sqrt(fs * fs + gs * gs)
        let cs = abs(fs) / d
        var r = copysign(d, f)
        let sn = gs / r
        r *= u

        return (cs, sn, r)
    }
}