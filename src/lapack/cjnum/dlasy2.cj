package cjnum.lapack.cjnum

import cjnum.blas.blas64.*
import std.math.*

extend LFloat64Implementation {
    public func dlasy2(tranl: Bool, tranr: Bool, isgn: Int64, n1: Int64, n2: Int64, tl: Array<Float64>, ldtl: Int64,
        tr: Array<Float64>, ldtr: Int64, b: Array<Float64>, ldb: Int64, x: Array<Float64>, ldx: Int64): (Float64, Float64, 
        Bool) {
        var scale: Float64 = 0.0
        var xnorm: Float64 = 0.0
        var ok: Bool = true

        if (n1 == 0 || n2 == 0) {
            return (scale, xnorm, ok)
        }

        let eps = dlamchP
        let smlnum = dlamchS / eps
        let sgn = Float64(isgn)

        if (n1 == 1 && n2 == 1) {
            var tau1 = tl[0] + sgn * tr[0]
            var bet = abs(tau1)
            if (bet <= smlnum) {
                tau1 = smlnum
                bet = smlnum
                ok = false
            }
            scale = 1.0
            let gam = abs(b[0])
            if (smlnum * gam > bet) {
                scale = 1.0 / gam
            }
            x[0] = b[0] * scale / tau1
            xnorm = abs(x[0])
            return (scale, xnorm, ok)
        }

        if (n1 + n2 == 3) {
            var smin: Float64
            var tmp: Array<Float64> = Array<Float64>(4, item: 0.0)
            var btmp: Array<Float64> = Array<Float64>(2, item: 0.0)
            if (n1 == 1 && n2 == 2) {
                smin = abs(tl[0])
                smin = max(smin, max(abs(tr[0]), abs(tr[1])))
                smin = max(smin, max(abs(tr[ldtr]), abs(tr[ldtr + 1])))
                smin = max(eps * smin, smlnum)
                tmp[0] = tl[0] + sgn * tr[0]
                tmp[3] = tl[0] + sgn * tr[ldtr + 1]
                if (tranr) {
                    tmp[1] = sgn * tr[1]
                    tmp[2] = sgn * tr[ldtr]
                } else {
                    tmp[1] = sgn * tr[ldtr]
                    tmp[2] = sgn * tr[1]
                }
                btmp[0] = b[0]
                btmp[1] = b[1]
            } else {
                smin = abs(tr[0])
                smin = max(smin, max(abs(tl[0]), abs(tl[1])))
                smin = max(smin, max(abs(tl[ldtl]), abs(tl[ldtl + 1])))
                smin = max(eps * smin, smlnum)
                tmp[0] = tl[0] + sgn * tr[0]
                tmp[3] = tl[ldtl + 1] + sgn * tr[0]
                if (tranl) {
                    tmp[1] = tl[ldtl]
                    tmp[2] = tl[1]
                } else {
                    tmp[1] = tl[1]
                    tmp[2] = tl[ldtl]
                }
                btmp[0] = b[0]
                btmp[1] = b[ldb]
            }

            let nFloat64 = nFloat64Implementation()
            let ipiv = nFloat64.idamax(tmp.size, tmp[..], 1)

            var u11 = tmp[ipiv]
            if (abs(u11) <= smin) {
                ok = false
                u11 = smin
            }
            let locu12 = Array<Int64>([1, 0, 3, 2])
            let u12 = tmp[locu12[ipiv]]
            let locl21 = Array<Int64>([2, 3, 0, 1])
            let l21 = tmp[locl21[ipiv]] / u11
            let locu22 = Array<Int64>([3, 2, 1, 0])
            var u22 = tmp[locu22[ipiv]] - l21 * u12
            if (abs(u22) <= smin) {
                ok = false
                u22 = smin
            }
            if ((ipiv & 0x2) != 0) {
                (btmp[0], btmp[1]) = (btmp[1], btmp[0] - l21 * btmp[1])
            } else {
                btmp[1] -= l21 * btmp[0]
            }
            scale = 1.0
            if (2.0 * smlnum * abs(btmp[1]) > abs(u22) || 2.0 * smlnum * abs(btmp[0]) > abs(u11)) {
                scale = 0.5 / max(abs(btmp[0]), abs(btmp[1]))
                btmp[0] *= scale
                btmp[1] *= scale
            }

            var x22 = btmp[1] / u22
            var x21 = btmp[0] / u11 - (u12 / u11) * x22
            if ((ipiv & 0x1) != 0) {
                (x21, x22) = (x22, x21)
            }
            x[0] = x21
            if (n1 == 1) {
                x[1] = x22
                xnorm = abs(x[0]) + abs(x[1])
            } else {
                x[ldx] = x22
                xnorm = max(abs(x[0]), abs(x[ldx]))
            }
            return (scale, xnorm, ok)
        }

        var smin = max(abs(tr[0]), abs(tr[1]))
        smin = max(smin, max(abs(tr[ldtr]), abs(tr[ldtr + 1])))
        smin = max(smin, max(abs(tl[0]), abs(tl[1])))
        smin = max(smin, max(abs(tl[ldtl]), abs(tl[ldtl + 1])))
        smin = max(eps * smin, smlnum)

        let t = Array<Array<Float64>>(
            [Array<Float64>(4, item: 0.0), Array<Float64>(4, item: 0.0), Array<Float64>(4, item: 0.0),
                Array<Float64>(4, item: 0.0)])
        t[0][0] = tl[0] + sgn * tr[0]
        t[1][1] = tl[0] + sgn * tr[ldtr + 1]
        t[2][2] = tl[ldtl + 1] + sgn * tr[0]
        t[3][3] = tl[ldtl + 1] + sgn * tr[ldtr + 1]
        if (tranl) {
            t[0][2] = tl[ldtl]
            t[1][3] = tl[ldtl]
            t[2][0] = tl[1]
            t[3][1] = tl[1]
        } else {
            t[0][2] = tl[1]
            t[1][3] = tl[1]
            t[2][0] = tl[ldtl]
            t[3][1] = tl[ldtl]
        }
        if (tranr) {
            t[0][1] = sgn * tr[1]
            t[1][0] = sgn * tr[ldtr]
            t[2][3] = sgn * tr[1]
            t[3][2] = sgn * tr[ldtr]
        } else {
            t[0][1] = sgn * tr[ldtr]
            t[1][0] = sgn * tr[1]
            t[2][3] = sgn * tr[ldtr]
            t[3][2] = sgn * tr[1]
        }

        let btmp = Array<Float64>(4, item: 0.0)
        btmp[0] = b[0]
        btmp[1] = b[1]
        btmp[2] = b[ldb]
        btmp[3] = b[ldb + 1]

        var jpiv = Array<Int64>(4, item: 0)
        for (i in 0..3) {
            var xmax: Float64 = 0.0
            var ipsv: Int64 = 0
            var jpsv: Int64 = 0
            for (ip in i..4) {
                for (jp in i..4) {
                    if (abs(t[ip][jp]) >= xmax) {
                        xmax = abs(t[ip][jp])
                        ipsv = ip
                        jpsv = jp
                    }
                }
            }
            if (ipsv != i) {
                (t[ipsv], t[i]) = (t[i], t[ipsv])
                (btmp[ipsv], btmp[i]) = (btmp[i], btmp[ipsv])
            }
            if (jpsv != i) {
                for (k in 0..4) {
                    (t[k][jpsv], t[k][i]) = (t[k][i], t[k][jpsv])
                }
            }
            jpiv[i] = jpsv
            if (abs(t[i][i]) < smin) {
                ok = false
                t[i][i] = smin
            }
            for (k in i + 1..4) {
                t[k][i] /= t[i][i]
                btmp[k] -= t[k][i] * btmp[i]
                for (j in i + 1..4) {
                    t[k][j] -= t[k][i] * t[i][j]
                }
            }
        }
        if (abs(t[3][3]) < smin) {
            ok = false
            t[3][3] = smin
        }
        scale = 1.0
        if (8.0 * smlnum * abs(btmp[0]) > abs(t[0][0]) || 8.0 * smlnum * abs(btmp[1]) > abs(t[1][1]) ||
            8.0 * smlnum * abs(btmp[2]) > abs(t[2][2]) || 8.0 * smlnum * abs(btmp[3]) > abs(t[3][3])) {
            var maxbtmp = max(abs(btmp[0]), abs(btmp[1]))
            maxbtmp = max(maxbtmp, max(abs(btmp[2]), abs(btmp[3])))
            scale = (1.0 / 8.0) / maxbtmp
            btmp[0] *= scale
            btmp[1] *= scale
            btmp[2] *= scale
            btmp[3] *= scale
        }

        var tmp: Array<Float64> = Array<Float64>(4, item: 0.0)
        for (i in 3..=0 : -1) {
            let temp = 1.0 / t[i][i]
            tmp[i] = btmp[i] * temp
            for (j in i + 1..4) {
                tmp[i] -= temp * t[i][j] * tmp[j]
            }
        }
        for (i in 2..=0 : -1) {
            if (jpiv[i] != i) {
                (tmp[i], tmp[jpiv[i]]) = (tmp[jpiv[i]], tmp[i])
            }
        }
        x[0] = tmp[0]
        x[1] = tmp[1]
        x[ldx] = tmp[2]
        x[ldx + 1] = tmp[3]
        xnorm = max(abs(tmp[0]) + abs(tmp[1]), abs(tmp[2]) + abs(tmp[3]))
        return (scale, xnorm, ok)
    }
}
