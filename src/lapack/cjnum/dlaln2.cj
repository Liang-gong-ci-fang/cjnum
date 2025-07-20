package cjnum.lapack.cjnum

import cjnum.complex.*
import std.math.*

extend LFloat64Implementation {
    public func dlaln2(trans: Bool, na: Int64, nw: Int64, smin: Float64, ca: Float64, a: Array<Float64>, lda: Int64,
        d1: Float64, d2: Float64, b: Array<Float64>, ldb: Int64, wr: Float64, wi: Float64, x: Array<Float64>, ldx: Int64): (Float64, 
        Float64, Bool) {
        var scale: Float64 = 0.0
        var xnorm: Float64 = 0.0
        var ok: Bool = false

        if (na != 1 && na != 2) {
            throw Exception(badNa)
        }
        if (nw != 1 && nw != 2) {
            throw Exception(badNw)
        }
        if (lda < na) {
            throw Exception(badLdA)
        }
        if (a.size < (na - 1) * lda + na) {
            throw Exception(shortA)
        }
        if (ldb < nw) {
            throw Exception(badLdB)
        }
        if (b.size < (na - 1) * ldb + nw) {
            throw Exception(shortB)
        }
        if (ldx < nw) {
            throw Exception(badLdX)
        }
        if (x.size < (na - 1) * ldx + nw) {
            throw Exception(shortX)
        }

        let smlnum = 2.0 * dlamchS
        let bignum = 1.0 / smlnum
        let smini = max(smin, smlnum)

        ok = true
        scale = 1.0

        if (na == 1) {
            if (nw == 1) {
                var csr = ca * a[0] - wr * d1
                var cnorm = abs(csr)

                if (cnorm < smini) {
                    csr = smini
                    cnorm = smini
                    ok = false
                }

                let bnorm = abs(b[0])
                if (cnorm < 1.0 && bnorm > max(1.0, bignum * cnorm)) {
                    scale = 1.0 / bnorm
                }

                x[0] = b[0] * scale / csr
                xnorm = abs(x[0])

                return (scale, xnorm, ok)
            }

            var csr = ca * a[0] - wr * d1
            var csi = -wi * d1
            var cnorm = abs(csr) + abs(csi)

            if (cnorm < smini) {
                csr = smini
                csi = 0.0
                cnorm = smini
                ok = false
            }

            let bnorm = abs(b[0]) + abs(b[1])
            if (cnorm < 1.0 && bnorm > max(1.0, bignum * cnorm)) {
                scale = 1.0 / bnorm
            }

            let cx = Complex128(scale * b[0], scale * b[1]) / Complex128(csr, csi)
            (x[0], x[1]) = (cx.real, cx.imag)
            xnorm = abs(x[0]) + abs(x[1])

            return (scale, xnorm, ok)
        }

        let crv = [
                ca * a[0] - wr * d1,
                ca * a[1],
                ca * a[lda],
                ca * a[lda + 1] - wr * d2
            ]
        
        if (trans) {
            crv[1] = ca * a[lda]
            crv[2] = ca * a[1]
        }

        let pivot = [
                [0, 1, 2, 3],
                [1, 0, 3, 2],
                [2, 3, 0, 1],
                [3, 2, 1, 0]
            ]
        
        if (nw == 1) {
            var cmax = 0.0
            var icmax = 0
            for (j in 0..crv.size) {
                let v = abs(crv[j])
                if (v > cmax) {
                    cmax = v
                    icmax = j
                }
            }

            if (cmax < smini) {
                let bnorm = max(abs(b[0]), abs(b[ldb]))
                if (smini < 1.0 && bnorm > max(1.0, bignum * smini)) {
                    scale = 1.0 / bnorm
                }
                let temp = scale / smini
                x[0] = temp * b[0]
                x[ldx] = temp * b[ldb]
                xnorm = temp * bnorm
                ok = false

                return (scale, xnorm, ok)
            }

            let ur11 = crv[icmax]
            let ur12 = crv[pivot[icmax][1]]
            let cr21 = crv[pivot[icmax][2]]
            let cr22 = crv[pivot[icmax][3]]
            let ur11r = 1.0 / ur11
            let lr21 = ur11r * cr21
            var ur22 = cr22 - ur12 * lr21

            if (abs(ur22) < smini) {
                ur22 = smini
                ok = false
            }

            var br1: Float64
            var br2: Float64
            if (icmax > 1) {
                br1 = b[ldb]
                br2 = b[0]
            } else {
                br1 = b[0]
                br2 = b[ldb]
            }
            br2 -= lr21 * br1

            let bbnd = max(abs(ur22 * ur11r * br1), abs(br2))
            if (bbnd > 1.0 && abs(ur22) < 1.0 && bbnd >= bignum * abs(ur22)) {
                scale = 1.0 / bbnd
            }

            let xr2 = br2 * scale / ur22
            let xr1 = scale * br1 * ur11r - ur11r * ur12 * xr2
            if (icmax % 2 != 0) {
                x[0] = xr2
                x[ldx] = xr1
            } else {
                x[0] = xr1
                x[ldx] = xr2
            }
            xnorm = max(abs(xr1), abs(xr2))

            if (xnorm > 1.0 && cmax > 1.0 && xnorm > bignum / cmax) {
                let temp = cmax / bignum
                x[0] *= temp
                x[ldx] *= temp
                xnorm *= temp
                scale *= temp
            }

            return (scale, xnorm, ok)
        }

        let civ = [
                -wi * d1,
                0.0,
                0.0,
                -wi * d2
            ]
        var cmax = 0.0
        var icmax = 0
        for (j in 0..crv.size) {
            let v = abs(crv[j])
            if (v + abs(civ[j]) > cmax) {
                cmax = v + abs(civ[j])
                icmax = j
            }
        }

        if (cmax < smini) {
            let br1 = abs(b[0]) + abs(b[1])
            let br2 = abs(b[ldb]) + abs(b[ldb + 1])
            let bnorm = max(br1, br2)
            if (smini < 1.0 && bnorm > 1.0 && bnorm > bignum * smini) {
                scale = 1.0 / bnorm
            }
            let temp = scale / smini
            x[0] = temp * b[0]
            x[1] = temp * b[1]
            x[ldb] = temp * b[ldb]
            x[ldb + 1] = temp * b[ldb + 1]
            xnorm = temp * bnorm
            ok = false

            return (scale, xnorm, ok)
        }

        let ur11 = crv[icmax]
        let ui11 = civ[icmax]
        let ur12 = crv[pivot[icmax][1]]
        let ui12 = civ[pivot[icmax][1]]
        let cr21 = crv[pivot[icmax][2]]
        let ci21 = civ[pivot[icmax][2]]
        let cr22 = crv[pivot[icmax][3]]
        let ci22 = civ[pivot[icmax][3]]
        var ur11r: Float64
        var ui11r: Float64
        var lr21: Float64
        var li21: Float64
        var ur12s: Float64
        var ui12s: Float64
        var ur22: Float64
        var ui22: Float64
        if (icmax == 0 || icmax == 3) {
            if (abs(ur11) > abs(ui11)) {
                let temp = ui11 / ur11
                ur11r = 1.0 / (ur11 * (1.0 + temp * temp))
                ui11r = -temp * ur11r
            } else {
                let temp = ur11 / ui11
                ui11r = -1.0 / (ui11 * (1.0 + temp * temp))
                ur11r = -temp * ui11r
            }
            lr21 = cr21 * ur11r
            li21 = cr21 * ui11r
            ur12s = ur12 * ur11r
            ui12s = ur12 * ui11r
            ur22 = cr22 - ur12 * lr21
            ui22 = ci22 - ur12 * li21
        } else {
            ur11r = 1.0 / ur11
            ui11r = 0.0
            lr21 = cr21 * ur11r
            li21 = ci21 * ur11r
            ur12s = ur12 * ur11r
            ui12s = ui12 * ur11r
            ur22 = cr22 - ur12 * lr21 + ui12 * li21
            ui22 = -ur12 * li21 - ui12 * lr21
        }
        let u22abs = abs(ur22) + abs(ui22)

        if (u22abs < smini) {
            ur22 = smini
            ui22 = 0.0
            ok = false
        }

        var br1: Float64
        var bi1: Float64
        var br2: Float64
        var bi2: Float64
        if (icmax > 1) {
            br1 = b[ldb]
            bi1 = b[ldb + 1]
            br2 = b[0]
            bi2 = b[1]
        } else {
            br1 = b[0]
            bi1 = b[1]
            br2 = b[ldb]
            bi2 = b[ldb + 1]
        }
        br2 += -lr21 * br1 + li21 * bi1
        bi2 += -li21 * br1 - lr21 * bi1

        let bbnd1 = u22abs * (abs(ur11r) + abs(ui11r)) * (abs(br1) + abs(bi1))
        let bbnd2 = abs(br2) + abs(bi2)
        let bbnd = max(bbnd1, bbnd2)
        if (bbnd > 1.0 && u22abs < 1.0 && bbnd >= bignum * u22abs) {
            scale = 1.0 / bbnd
            br1 *= scale
            bi1 *= scale
            br2 *= scale
            bi2 *= scale
        }

        let cx2 = Complex128(br2, bi2) / Complex128(ur22, ui22)
        let (xr2, xi2) = (cx2.real, cx2.imag)
        let xr1 = ur11r * br1 - ui11r * bi1 - ur12s * xr2 + ui12s * xi2
        let xi1 = ui11r * br1 + ur11r * bi1 - ui12s * xr2 - ur12s * xi2
        if (icmax % 2 != 0) {
            x[0] = xr2
            x[1] = xi2
            x[ldx] = xr1
            x[ldx + 1] = xi1
        } else {
            x[0] = xr1
            x[1] = xi1
            x[ldx] = xr2
            x[ldx + 1] = xi2
        }
        xnorm = max(abs(xr1) + abs(xi1), abs(xr2) + abs(xi2))

        if (xnorm > 1.0 && cmax > 1.0 && xnorm > bignum / cmax) {
            let temp = cmax / bignum
            x[0] *= temp
            x[1] *= temp
            x[ldx] *= temp
            x[ldx + 1] *= temp
            xnorm *= temp
            scale *= temp
        }

        return (scale, xnorm, ok)
    }
}
