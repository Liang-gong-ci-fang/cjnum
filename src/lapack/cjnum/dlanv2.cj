package cjnum.lapack.cjnum

import std.math.*

extend LFloat64Implementation {
    public func dlanv2(a: Float64, b: Float64, c: Float64, d: Float64): (Float64, Float64, Float64, Float64, Float64, 
        Float64, Float64, Float64, Float64, Float64) {
        var aa: Float64
        var bb: Float64
        var cc: Float64
        var dd: Float64
        var cs: Float64
        var sn: Float64
        var rt1r: Float64
        var rt2r: Float64
        var rt1i: Float64 = 0.0
        var rt2i: Float64 = 0.0

        if (c == 0.0) {
            aa = a
            bb = b
            cc = 0.0
            dd = d
            cs = 1.0
            sn = 0.0
        } else if (b == 0.0) {
            aa = d
            bb = -c
            cc = 0.0
            dd = a
            cs = 0.0
            sn = 1.0
        } else if (a == d && b * c < 0.0) {
            aa = a
            bb = b
            cc = c
            dd = d
            cs = 1.0
            sn = 0.0
        } else {
            var temp = a - d
            var p = temp / 2.0
            var bcmax = max(abs(b), abs(c))
            var bcmis = min(abs(b), abs(c))
            if (b * c < 0.0) {
                bcmis *= -1.0
            }
            var scale = max(abs(p), bcmax)
            var z = p / scale * p + bcmax / scale * bcmis
            let eps = dlamchP

            if (z >= 4.0 * eps) {
                if (p > 0.0) {
                    z = p + sqrt(scale) * sqrt(z)
                } else {
                    z = p - sqrt(scale) * sqrt(z)
                }
                aa = d + z
                dd = d - bcmax / z * bcmis
                let tau = dlapy2(c, z)
                cs = z / tau
                sn = c / tau
                bb = b - c
                cc = 0.0
            } else {
                let safmn2 = pow(Float64(dlamchB), log(dlamchS / dlamchE) / log(Float64(dlamchB)) / 2.0)
                let safmx2 = 1.0 / safmn2
                var sigma = b + c

                for (iter in 0..20) {
                    scale = max(abs(temp), abs(sigma))
                    if (scale >= safmx2) {
                        sigma *= safmn2
                        temp *= safmn2
                    } else if (scale <= safmn2) {
                        sigma *= safmx2
                        temp *= safmx2
                    } else {
                        break
                    }
                }
                p = temp / 2.0
                var tau = dlapy2(sigma, temp)
                cs = sqrt((1.0 + abs(sigma) / tau) / 2.0)
                sn = -p / (tau * cs)
                if (sigma < 0.0) {
                    sn *= -1.0
                }

                aa = a * cs + b * sn
                bb = -a * sn + b * cs
                cc = c * cs + d * sn
                dd = -c * sn + d * cs

                let a_ = aa * cs + cc * sn
                let b_ = bb * cs + dd * sn
                let c_ = -aa * sn + cc * cs
                let d_ = -bb * sn + dd * cs

                temp = (a_ + d_) / 2.0
                aa = temp
                bb = b_
                cc = c_
                dd = temp

                if (cc != 0.0) {
                    if (bb != 0.0) {
                        if (bb * cc > 0.0) {
                            let sab = sqrt(abs(bb))
                            let sac = sqrt(abs(cc))
                            p = sab * sac
                            if (cc < 0.0) {
                                p *= -1.0
                            }
                            tau = 1.0 / sqrt(abs(bb + cc))
                            aa = temp + p
                            bb = bb - cc
                            cc = 0.0
                            dd = temp - p
                            let cs1 = sab * tau
                            let sn1 = sac * tau
                            (cs, sn) = (cs * cs1 - sn * sn1, cs * sn1 + sn * cs1)
                        }
                    } else {
                        bb = -cc
                        cc = 0.0
                        (cs, sn) = (-sn, cs)
                    }
                }
            }
        }

        rt1r = aa
        rt2r = dd
        if (cc != 0.0) {
            rt1i = sqrt(abs(bb)) * sqrt(abs(cc))
            rt2i = -rt1i
        }
        return (aa, bb, cc, dd, rt1r, rt1i, rt2r, rt2i, cs, sn)
    }
}
