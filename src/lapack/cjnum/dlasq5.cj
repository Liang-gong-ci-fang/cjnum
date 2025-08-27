package cjnum.lapack.cjnum

extend LFloat64Implementation {
    public func dlasq5(i0: Int64, n0: Int64, z: Array<Float64>, pp: Int64, tau: Float64, sigma: Float64): (Int64, Int64, 
        Int64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64) {
        var tau_ = tau

        var dmin: Float64 = 0.0
        var dmin1: Float64 = 0.0
        var dmin2: Float64 = 0.0
        var dn: Float64 = 0.0
        var dnm1: Float64 = 0.0
        var dnm2: Float64 = 0.0

        if (i0 < 0) {
            throw Exception(i0LT0)
        } else if (n0 < 0) {
            throw Exception(n0LT0)
        } else if (z.size < 4 * n0) {
            throw Exception(shortZ)
        } else if (pp != 0 && pp != 1) {
            throw Exception(badPp)
        }

        if (n0 - i0 - 1 <= 0) {
            return (i0, n0, pp, tau_, sigma, dmin, dmin1, dmin2, dn, dnm1, dnm2)
        }

        let eps = dlamchP
        let dthresh = eps * (sigma + tau_)
        if (tau_ < dthresh * 0.5) {
            tau_ = 0.0
        }

        var j4 = 0
        var emin = 0.0

        if (tau_ != 0.0) {
            j4 = 4 * i0 + pp
            emin = z[j4 + 4]
            var d = z[j4] - tau_
            dmin = d

            if (pp == 0) {
                for (j4loop in 4 * (i0 + 1)..=4 * ((n0 + 1) - 3) : 4) {
                    let j4 = j4loop - 1
                    z[j4 - 2] = d + z[j4 - 1]
                    let tmp = z[j4 + 1] / z[j4 - 2]
                    d = d * tmp - tau_
                    dmin = min(dmin, d)
                    z[j4] = z[j4 - 1] * tmp
                    emin = min(z[j4], emin)
                }
            } else {
                for (j4loop in 4 * (i0 + 1)..=4 * ((n0 + 1) - 3) : 4) {
                    let j4 = j4loop - 1
                    z[j4 - 3] = d + z[j4]
                    let tmp = z[j4 + 2] / z[j4 - 3]
                    d = d * tmp - tau_
                    dmin = min(dmin, d)
                    z[j4 - 1] = z[j4] * tmp
                    emin = min(z[j4 - 1], emin)
                }
            }

            dnm2 = d
            dmin2 = dmin
            j4 = 4 * ((n0 + 1) - 2) - pp - 1
            let j4p2 = j4 + 2 * pp - 1
            z[j4 - 2] = dnm2 + z[j4p2]
            z[j4] = z[j4p2 + 2] * (z[j4p2] / z[j4 - 2])
            dnm1 = z[j4p2 + 2] * (dnm2 / z[j4 - 2]) - tau_
            dmin = min(dmin, dnm1)

            dmin1 = dmin
            j4 += 4
            let j4p2_next = j4 + 2 * pp - 1
            z[j4 - 2] = dnm1 + z[j4p2_next]
            z[j4] = z[j4p2_next + 2] * (z[j4p2_next] / z[j4 - 2])
            dn = z[j4p2_next + 2] * (dnm1 / z[j4 - 2]) - tau_
            dmin = min(dmin, dn)
        } else {
            j4 = 4 * (i0 + 1) + pp - 4
            emin = z[j4 + 4]
            var d = z[j4] - tau_
            dmin = d

            if (pp == 0) {
                for (j4loop in 4 * (i0 + 1)..=4 * ((n0 + 1) - 3) : 4) {
                    let j4 = j4loop - 1
                    z[j4 - 2] = d + z[j4 - 1]
                    let tmp = z[j4 + 1] / z[j4 - 2]
                    d = d * tmp - tau_
                    if (d < dthresh) {
                        d = 0.0
                    }
                    dmin = min(dmin, d)
                    z[j4] = z[j4 - 1] * tmp
                    emin = min(z[j4], emin)
                }
            } else {
                for (j4loop in 4 * (i0 + 1)..=4 * ((n0 + 1) - 3) : 4) {
                    let j4 = j4loop - 1
                    z[j4 - 3] = d + z[j4]
                    let tmp = z[j4 + 2] / z[j4 - 3]
                    d = d * tmp - tau_
                    if (d < dthresh) {
                        d = 0.0
                    }
                    dmin = min(dmin, d)
                    z[j4 - 1] = z[j4] * tmp
                    emin = min(z[j4 - 1], emin)
                }
            }

            dnm2 = d
            dmin2 = dmin
            j4 = 4 * ((n0 + 1) - 2) - pp - 1
            let j4p2 = j4 + 2 * pp - 1
            z[j4 - 2] = dnm2 + z[j4p2]
            z[j4] = z[j4p2 + 2] * (z[j4p2] / z[j4 - 2])
            dnm1 = z[j4p2 + 2] * (dnm2 / z[j4 - 2]) - tau_
            dmin = min(dmin, dnm1)

            dmin1 = dmin
            j4 += 4
            let j4p2_next = j4 + 2 * pp - 1
            z[j4 - 2] = dnm1 + z[j4p2_next]
            z[j4] = z[j4p2_next + 2] * (z[j4p2_next] / z[j4 - 2])
            dn = z[j4p2_next + 2] * (dnm1 / z[j4 - 2]) - tau_
            dmin = min(dmin, dn)
        }
        z[j4 + 2] = dn
        z[4 * (n0 + 1) - pp - 1] = emin
        return (i0, n0, pp, tau_, sigma, dmin, dmin1, dmin2, dn, dnm1, dnm2)
    }
}
