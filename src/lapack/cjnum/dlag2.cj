package cjnum.lapack.cjnum

import std.math.*
import cjnum.util.*

extend LFloat64Implementation {
    public func dlag2(a: Array<Float64>, lda: Int64, b: Array<Float64>, ldb: Int64): (Float64, Float64, Float64, 
        Float64, Float64) {
        var scale1: Float64 = 0.0
        var scale2: Float64 = 0.0
        var wr1: Float64 = 0.0
        var wr2: Float64 = 0.0
        var wi: Float64 = 0.0
        if (lda < 2) {
            throw Exception(badLdA)
        } else if (ldb < 2) {
            throw Exception(badLdB)
        } else if (a.size < lda + 2) {
            throw Exception(shortA)
        } else if (b.size < ldb + 2) {
            throw Exception(shortB)
        }

        const safmin = dlamchS
        const safmax = 1.0 / safmin
        const fuzzy1 = 1.0 + 1e-5

        let rtmin = sqrt(safmin)
        let rtmax = 1.0 / rtmin

        var anorm = nanMax(abs(a[0]) + abs(a[lda]), abs(a[1]) + abs(a[lda + 1]))
        anorm = nanMax(anorm, safmin)
        let ascale = 1.0 / anorm
        let a11 = ascale * a[0]
        let a21 = ascale * a[lda]
        let a12 = ascale * a[1]
        let a22 = ascale * a[lda + 1]

        var b11 = b[0]
        var b12 = b[1]
        var b22 = b[ldb + 1]
        let bmin = rtmin * nanMax(nanMax(abs(b11), abs(b12)), nanMax(abs(b22), rtmin))
        if (abs(b11) < bmin) {
            b11 = copysign(bmin, b11)
        }
        if (abs(b22) < bmin) {
            b22 = copysign(bmin, b22)
        }

        let bnorm = nanMax(nanMax(abs(b11), abs(b12) + abs(b22)), safmin)
        let bsize = nanMax(abs(b11), abs(b22))
        let bscale = 1.0 / bsize
        b11 *= bscale
        b12 *= bscale
        b22 *= bscale

        var as12: Float64
        var abi22: Float64
        var pp: Float64
        var qq: Float64
        var shift: Float64
        let binv11 = 1.0 / b11
        let binv22 = 1.0 / b22
        let s1 = a11 * binv11
        let s2 = a22 * binv22
        if (abs(s1) <= abs(s2)) {
            shift = s1
            as12 = a12 - shift * b12
            let as22 = a22 - shift * b22
            let ss = a21 * (binv11 * binv22)
            abi22 = as22 * binv22 - ss * b12
            pp = 0.5 * abi22
            qq = ss * as12
        } else {
            shift = s2
            as12 = a12 - shift * b12
            let as11 = a11 - shift * b11
            let ss = a21 * (binv11 * binv22)
            abi22 = -ss * b12
            pp = 0.5 * (as11 * binv11 + abi22)
            qq = ss * as12
        }
        var discr: Float64
        var r: Float64
        if (abs(pp * rtmin) >= 1.0) {
            let tmp = rtmin * pp
            discr = tmp * tmp + qq * safmin
            r = sqrt(abs(discr)) * rtmax
        } else {
            let pp2 = pp * pp
            if (pp2 + abs(qq) <= safmin) {
                let tmp = rtmax * pp
                discr = tmp * tmp + qq * safmax
                r = sqrt(abs(discr)) * rtmin
            } else {
                discr = pp2 + qq
                r = sqrt(abs(discr))
            }
        }

        if (discr >= 0.0 || r == 0.0) {
            let sum = pp + copysign(r, pp)
            let diff = pp - copysign(r, pp)
            let wbig = shift + sum

            var wsmall = shift + diff
            if (0.5 * abs(wbig) > nanMax(abs(wsmall), safmin)) {
                let wdet = (a11 * a22 - a12 * a21) * (binv11 * binv22)
                wsmall = wdet / wbig
            }
            if (pp > abi22) {
                wr1 = nanMin(wbig, wsmall)
                wr2 = nanMax(wbig, wsmall)
            } else {
                wr1 = nanMax(wbig, wsmall)
                wr2 = nanMin(wbig, wsmall)
            }
        } else {
            wr1 = shift + pp
            wr2 = wr1
            wi = r
        }

        let c1 = bsize * (safmin * nanMax(1.0, ascale))
        let c2 = safmin * nanMax(1.0, bnorm)
        let c3 = bsize * safmin
        var c4 = 1.0
        var c5 = 1.0
        if (ascale <= 1.0 || bsize <= 1.0) {
            c5 = nanMin(1.0, ascale * bsize)
            if (ascale <= 1.0 && bsize <= 1.0) {
                c4 = nanMin(1.0, (ascale / safmin) * bsize)
            }
        }

        let wabs = abs(wr1) + abs(wi)
        var wsize = nanMax(nanMax(safmin, c1), nanMax(fuzzy1 * (wabs * c2 + c3), nanMin(c4, 0.5 * nanMax(wabs, c5))))
        let maxABsize = nanMax(ascale, bsize)
        let minABsize = nanMin(ascale, bsize)
        if (wsize != 1.0) {
            let wscale = 1.0 / wsize
            if (wsize > 1.0) {
                scale1 = (maxABsize * wscale) * minABsize
            } else {
                scale1 = (minABsize * wscale) * maxABsize
            }
            wr1 *= wscale
            if (wi != 0.0) {
                wi *= wscale
                wr2 = wr1
                scale2 = scale1
            }
        } else {
            scale1 = ascale * bsize
            scale2 = scale1
        }

        if (wi == 0.0) {
            wsize = nanMax(nanMax(safmin, c1), nanMax(fuzzy1 * (abs(wr2) * c2 + c3), nanMin(c4, 0.5 * nanMax(abs(wr2), c5))))
            if (wsize != 1.0) {
                let wscale = 1.0 / wsize
                if (wsize > 1.0) {
                    scale2 = (maxABsize * wscale) * minABsize
                } else {
                    scale2 = (minABsize * wscale) * maxABsize
                }
                wr2 *= wscale
            } else {
                scale2 = ascale * bsize
            }
        }

        return (scale1, scale2, wr1, wr2, wi)
    }
}
