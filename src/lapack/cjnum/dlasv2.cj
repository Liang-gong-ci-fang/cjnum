package cjnum.lapack.cjnum

import std.math.*
import cjnum.util.*

extend LFloat64Implementation {
    public func dlasv2(f: Float64, g: Float64, h: Float64): (ssmin: Float64, ssmax: Float64, snr: Float64, csr: Float64, 
        snl: Float64, csl: Float64) {
        var ssmin: Float64 = 0.0
        var ssmax: Float64 = 0.0
        var snr: Float64 = 0.0
        var csr: Float64 = 0.0
        var snl: Float64 = 0.0
        var csl: Float64 = 0.0

        var ft = f
        var fa = abs(ft)
        var ht = h
        var ha = abs(h)

        var pmax = 1
        let swap = ha > fa
        if (swap) {
            pmax = 3
            (ft, ht) = (ht, ft)
            (fa, ha) = (ha, fa)
        }
        let gt = g
        let ga = abs(gt)
        var clt: Float64 = 0.0
        var crt: Float64 = 0.0
        var slt: Float64 = 0.0
        var srt: Float64 = 0.0
        if (ga == 0.0) {
            ssmin = ha
            ssmax = fa
            clt = 1.0
            crt = 1.0
            slt = 0.0
            srt = 0.0
        } else {
            var gasmall = true
            if (ga > fa) {
                pmax = 2
                if ((fa / ga) < dlamchE) {
                    gasmall = false
                    ssmax = ga
                    if (ha > 1.0) {
                        ssmin = fa / (ga / ha)
                    } else {
                        ssmin = (fa / ga) * ha
                    }
                    clt = 1.0
                    slt = ht / gt
                    srt = 1.0
                    crt = ft / gt
                }
            }
            if (gasmall) {
                let d = fa - ha
                var l = d / fa
                if (d == fa) {
                    l = 1.0
                }
                let m = gt / ft
                var t = 2.0 - l
                let s = hypot(t, m)
                var r: Float64
                if (l == 0.0) {
                    r = abs(m)
                } else {
                    r = hypot(l, m)
                }
                let a = 0.5 * (s + r)
                ssmin = ha / a
                ssmax = fa * a
                if (m == 0.0) {
                    if (l == 0.0) {
                        t = copysign(2.0, ft) * copysign(1.0, gt)
                    } else {
                        t = gt / copysign(d, ft) + m / t
                    }
                } else {
                    t = (m / (s + t) + m / (r + l)) * (1.0 + a)
                }
                l = hypot(t, 2.0)
                crt = 2.0 / l
                srt = t / l
                clt = (crt + srt * m) / a
                slt = (ht / ft) * srt / a
            }
        }
        if (swap) {
            csl = srt
            snl = crt
            csr = slt
            snr = clt
        } else {
            csl = clt
            snl = slt
            csr = crt
            snr = srt
        }
        var tsign: Float64 = 0.0
        if (pmax == 1) {
            tsign = copysign(1.0, csr) * copysign(1.0, csl) * copysign(1.0, f)
        } else if (pmax == 2) {
            tsign = copysign(1.0, snr) * copysign(1.0, csl) * copysign(1.0, g)
        } else if (pmax == 3) {
            tsign = copysign(1.0, snr) * copysign(1.0, snl) * copysign(1.0, h)
        }
        ssmax = copysign(ssmax, tsign)
        ssmin = copysign(ssmin, tsign * copysign(1.0, f) * copysign(1.0, h))
        return (ssmin, ssmax, snr, csr, snl, csl)
    }
}
