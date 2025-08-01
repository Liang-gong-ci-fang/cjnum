package cjnum.lapack.cjnum

import cjnum.util.*
import std.math.*

extend LFloat64Implementation {
    public func dlasq4(i0: Int64, n0: Int64, z: Array<Float64>, pp: Int64, n0in: Int64, dmin: Float64, dmin1: Float64,
        dmin2: Float64, dn: Float64, dn1: Float64, dn2: Float64, tau: Float64, ttype: Int64, g: Float64): (Float64, 
        Int64, Float64) {
        var tau_ = tau
        var ttype_ = ttype
        var g_ = g
        if (i0 < 0) {
            throw Exception("i0LT0")
        }
        if (n0 < 0) {
            throw Exception("n0LT0")
        }
        if (z.size < 4 * n0) {
            throw Exception("shortZ")
        }
        if (pp != 0 && pp != 1) {
            throw Exception("badPp")
        }

        let cnst1 = 0.563
        let cnst2 = 1.01
        let cnst3 = 1.05
        let cnstthird = 0.333

        if (dmin <= 0.0) {
            tau_ = -dmin
            ttype_ = -1
            return (tau_, ttype_, g_)
        }
        var nn = 4 * (n0 + 1) + pp - 1
        var s = Float64.NaN
        if (n0in == n0) {
            if (dmin == dn || dmin == dn1) {
                var b1 = sqrt(z[nn - 3]) * sqrt(z[nn - 5])
                var b2 = sqrt(z[nn - 7]) * sqrt(z[nn - 9])
                var a2 = z[nn - 7] + z[nn - 5]
                if (dmin == dn && dmin1 == dn1) {
                    let gap2 = dmin2 - a2 - dmin2 / 4.0
                    var gap1: Float64
                    if (gap2 > 0.0 && gap2 > b2) {
                        gap1 = a2 - dn - (b2 / gap2) * b2
                    } else {
                        gap1 = a2 - dn - (b1 + b2)
                    }
                    if (gap1 > 0.0 && gap1 > b1) {
                        s = nanMax(dn - (b1 / gap1) * b1, 0.5 * dmin)
                        ttype_ = -2
                    } else {
                        s = 0.0
                        if (dn > b1) {
                            s = dn - b1
                        }
                        if (a2 > b1 + b2) {
                            s = min(s, a2 - (b1 + b2))
                        }
                        s = nanMax(s, cnstthird * dmin)
                        ttype_ = -3
                    }
                } else {
                    ttype_ = -4
                    s = dmin / 4.0
                    var gam: Float64
                    var np: Int64
                    if (dmin == dn) {
                        gam = dn
                        a2 = 0.0
                        if (z[nn - 5] > z[nn - 7]) {
                            return (tau_, ttype_, g_)
                        }
                        b2 = z[nn - 5] / z[nn - 7]
                        np = nn - 9
                    } else {
                        np = nn - 2 * pp
                        gam = dn1
                        if (z[np - 4] > z[np - 2]) {
                            return (tau_, ttype_, g_)
                        }
                        a2 = z[np - 4] / z[np - 2]
                        if (z[nn - 9] > z[nn - 11]) {
                            return (tau_, ttype_, g_)
                        }
                        b2 = z[nn - 9] / z[nn - 11]
                        np = nn - 13
                    }
                    a2 += b2
                    for (i4loop in np + 1..=4 * (i0 + 1) - 1 + pp : -4) {
                        let i4 = i4loop - 1
                        if (b2 == 0.0) {
                            break
                        }
                        b1 = b2
                        if (z[i4] > z[i4 - 2]) {
                            return (tau_, ttype_, g_)
                        }
                        b2 *= z[i4] / z[i4 - 2]
                        a2 += b2
                        if (100.0 * nanMax(b2, b1) < a2 || cnst1 < a2) {
                            break
                        }
                    }
                    a2 *= cnst3
                    if (a2 < cnst1) {
                        s = gam * (1.0 - sqrt(a2)) / (1.0 + a2)
                    }
                }
            } else if (dmin == dn2) {
                ttype_ = -5
                s = dmin / 4.0
                let np = nn - 2 * pp
                var b1 = z[np - 2]
                var b2 = z[np - 6]
                let gam = dn2
                if (z[np - 8] > b2 || z[np - 4] > b1) {
                    return (tau_, ttype_, g_)
                }
                var a2 = (z[np - 8] / b2) * (1.0 + z[np - 4] / b1)
                if (n0 - i0 > 2) {
                    b2 = z[nn - 13] / z[nn - 15]
                    a2 += b2
                    for (i4loop in (nn + 1) - 17..=4 * (i0 + 1) - 1 + pp : -4) {
                        let i4 = i4loop - 1
                        if (b2 == 0.0) {
                            break
                        }
                        b1 = b2
                        if (z[i4] > z[i4 - 2]) {
                            return (tau_, ttype_, g_)
                        }
                        b2 *= z[i4] / z[i4 - 2]
                        a2 += b2
                        if (100.0 * nanMax(b2, b1) < a2 || cnst1 < a2) {
                            break
                        }
                    }
                    a2 *= cnst3
                }
                if (a2 < cnst1) {
                    s = gam * (1.0 - sqrt(a2)) / (1.0 + a2)
                }
            } else {
                if (ttype_ == -6) {
                    g_ += cnstthird * (1.0 - g_)
                } else if (ttype_ == -18) {
                    g_ = cnstthird / 4.0
                } else {
                    g_ = 1.0 / 4.0
                }
                s = g_ * dmin
                ttype_ = -6
            }
        } else if (n0in == (n0 + 1)) {
            if (dmin1 == dn1 && dmin2 == dn2) {
                ttype_ = -7
                s = cnstthird * dmin1
                if (z[nn - 5] > z[nn - 7]) {
                    return (tau_, ttype_, g_)
                }
                var b1 = z[nn - 5] / z[nn - 7]
                var b2 = b1
                if (b2 != 0.0) {
                    for (i4loop in 4 * (n0 + 1) - 9 + pp..=4 * (i0 + 1) - 1 + pp : -4) {
                        let i4 = i4loop - 1
                        let a2 = b1
                        if (z[i4] > z[i4 - 2]) {
                            return (tau_, ttype_, g_)
                        }
                        b1 *= z[i4] / z[i4 - 2]
                        b2 += b1
                        if (100.0 * nanMax(b1, a2) < b2) {
                            break
                        }
                    }
                }
                b2 = sqrt(cnst3 * b2)
                let a2 = dmin1 / (1.0 + b2 * b2)
                let gap2 = 0.5 * dmin2 - a2
                if (gap2 > 0.0 && gap2 > b2 * a2) {
                    s = nanMax(s, a2 * (1.0 - cnst2 * a2 * (b2 / gap2) * b2))
                } else {
                    s = nanMax(s, a2 * (1.0 - cnst2 * b2))
                    ttype_ = -8
                }
            } else {
                s = dmin1 / 4.0
                if (dmin1 == dn1) {
                    s = 0.5 * dmin1
                }
                ttype_ = -9
            }
        } else if (n0in == (n0 + 2)) {
            if (dmin2 == dn2 && 2.0 * z[nn - 5] < z[nn - 7]) {
                ttype_ = -10
                s = cnstthird * dmin2
                if (z[nn - 5] > z[nn - 7]) {
                    return (tau_, ttype_, g_)
                }
                var b1 = z[nn - 5] / z[nn - 7]
                var b2 = b1
                if (b2 != 0.0) {
                    for (i4loop in 4 * (n0 + 1) - 9 + pp..=4 * (i0 + 1) - 1 + pp : -4) {
                        let i4 = i4loop - 1
                        if (z[i4] > z[i4 - 2]) {
                            return (tau_, ttype_, g_)
                        }
                        b1 *= z[i4] / z[i4 - 2]
                        b2 += b1
                        if (100.0 * b1 < b2) {
                            break
                        }
                    }
                }
                b2 = sqrt(cnst3 * b2)
                let a2 = dmin2 / (1.0 + b2 * b2)
                let gap2 = z[nn - 7] + z[nn - 9] - sqrt(z[nn - 11]) * sqrt(z[nn - 9]) - a2
                if (gap2 > 0.0 && gap2 > b2 * a2) {
                    s = nanMax(s, a2 * (1.0 - cnst2 * a2 * (b2 / gap2) * b2))
                } else {
                    s = nanMax(s, a2 * (1.0 - cnst2 * b2))
                }
            } else {
                s = dmin2 / 4.0
                ttype_ = -11
            }
        } else if (n0in > n0 + 2) {
            s = 0.0
            ttype_ = -12
        }
        tau_ = s
        return (tau_, ttype_, g_)
    }
}
