package cjnum.lapack.cjnum

extend LFloat64Implementation {
    func dlasq6(i0: Int64, n0: Int64, z: Array<Float64>, pp: Int64): (Float64, Float64, Float64, Float64, Float64, 
        Float64) {
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
            return (dmin, dmin1, dmin2, dn, dnm1, dnm2)
        }

        let safmin = dlamchS
        var j4 = 4 * (i0 + 1) + pp - 4
        var emin = z[j4 + 4]
        var d = z[j4]
        dmin = d

        if (pp == 0) {
            for (j4loop in 4 * (i0 + 1)..=4 * ((n0 + 1) - 3) : 4) {
                let j4 = j4loop - 1
                z[j4 - 2] = d + z[j4 - 1]
                if (z[j4 - 2] == 0.0) {
                    z[j4] = 0.0
                    d = z[j4 + 1]
                    dmin = d
                    emin = 0.0
                } else if (safmin * z[j4 + 1] < z[j4 - 2] && safmin * z[j4 - 2] < z[j4 + 1]) {
                    let tmp = z[j4 + 1] / z[j4 - 2]
                    z[j4] = z[j4 - 1] * tmp
                    d *= tmp
                } else {
                    z[j4] = z[j4 + 1] * (z[j4 - 1] / z[j4 - 2])
                    d = z[j4 + 1] * (d / z[j4 - 2])
                }
                dmin = min(dmin, d)
                emin = min(emin, z[j4])
            }
        } else {
            for (j4loop in 4 * (i0 + 1)..=4 * ((n0 + 1) - 3) : 4) {
                let j4 = j4loop - 1
                z[j4 - 3] = d + z[j4]
                if (z[j4 - 3] == 0.0) {
                    z[j4 - 1] = 0.0
                    d = z[j4 + 2]
                    dmin = d
                    emin = 0.0
                } else if (safmin * z[j4 + 2] < z[j4 - 3] && safmin * z[j4 - 3] < z[j4 + 2]) {
                    let tmp = z[j4 + 2] / z[j4 - 3]
                    z[j4 - 1] = z[j4] * tmp
                    d *= tmp
                } else {
                    z[j4 - 1] = z[j4 + 2] * (z[j4] / z[j4 - 3])
                    d = z[j4 + 2] * (d / z[j4 - 3])
                }
                dmin = min(dmin, d)
                emin = min(emin, z[j4 - 1])
            }
        }

        dnm2 = d
        dmin2 = dmin
        j4 = 4 * (n0 - 1) - pp - 1
        let j4p2 = j4 + 2 * pp - 1
        z[j4 - 2] = dnm2 + z[j4p2]
        if (z[j4 - 2] == 0.0) {
            z[j4] = 0.0
            dnm1 = z[j4p2 + 2]
            dmin = dnm1
            emin = 0.0
        } else if (safmin * z[j4p2 + 2] < z[j4 - 2] && safmin * z[j4 - 2] < z[j4p2 + 2]) {
            let tmp = z[j4p2 + 2] / z[j4 - 2]
            z[j4] = z[j4p2] * tmp
            dnm1 = dnm2 * tmp
        } else {
            z[j4] = z[j4p2 + 2] * (z[j4p2] / z[j4 - 2])
            dnm1 = z[j4p2 + 2] * (dnm2 / z[j4 - 2])
        }
        dmin = min(dmin, dnm1)
        dmin1 = dmin
        j4 += 4
        let j4p2_new = j4 + 2 * pp - 1
        z[j4 - 2] = dnm1 + z[j4p2_new]
        if (z[j4 - 2] == 0.0) {
            z[j4] = 0.0
            dn = z[j4p2_new + 2]
            dmin = dn
            emin = 0.0
        } else if (safmin * z[j4p2_new + 2] < z[j4 - 2] && safmin * z[j4 - 2] < z[j4p2_new + 2]) {
            let tmp = z[j4p2_new + 2] / z[j4 - 2]
            z[j4] = z[j4p2_new] * tmp
            dn = dnm1 * tmp
        } else {
            z[j4] = z[j4p2_new + 2] * (z[j4p2_new] / z[j4 - 2])
            dn = z[j4p2_new + 2] * (dnm1 / z[j4 - 2])
        }
        dmin = min(dmin, dn)
        z[j4 + 2] = dn
        z[4 * (n0 + 1) - pp - 1] = emin
        return (dmin, dmin1, dmin2, dn, dnm1, dnm2)
    }
}
