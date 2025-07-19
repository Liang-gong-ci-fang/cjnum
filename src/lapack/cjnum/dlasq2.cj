package cjnum.lapack.cjnum

import cjnum.lapack.*
import std.math.*

extend LFloat64Implementation {
    public func dlasq2(n: Int64, z: Array<Float64>): Int64 {
        var info: Int64 = 0
        if (n < 0) {
            throw Exception(nLT0)
        }

        if (n == 0) {
            return info
        }

        if (z.size < 4 * n) {
            throw Exception(shortZ)
        }

        if (n == 1) {
            if (z[0] < 0.0) {
                throw Exception(negZ)
            }
            return info
        }

        const cbias = 1.5

        let eps = dlamchP
        let safmin = dlamchS
        let tol = eps * 100.0
        let tol2 = tol * tol
        if (n == 2) {
            if (z[1] < 0.0 || z[2] < 0.0) {
                throw Exception(negZ)
            } else if (z[2] > z[0]) {
                (z[0], z[2]) = (z[2], z[0])
            }
            z[4] = z[0] + z[1] + z[2]
            if (z[1] > z[2] * tol2) {
                var t = 0.5 * (z[0] - z[2] + z[1])
                var s = z[2] * (z[1] / t)
                if (s <= t) {
                    s = z[2] * (z[1] / (t * (1.0 + sqrt(1.0 + s / t))))
                } else {
                    s = z[2] * (z[1] / (t + sqrt(t) * sqrt(t + s)))
                }
                t = z[0] + s + z[1]
                z[2] *= z[0] / t
                z[0] = t
            }
            z[1] = z[2]
            z[5] = z[1] + z[0]
            return info
        }
        z[2 * n - 1] = 0.0
        var emin = z[1]
        var d: Float64 = 0.0
        var e: Float64 = 0.0
        var qmax: Float64 = 0.0
        var i1: Int64
        var n1: Int64
        for (k in 0..2 * (n - 1) : 2) {
            if (z[k] < 0.0 || z[k + 1] < 0.0) {
                throw Exception(negZ)
            }
            d += z[k]
            e += z[k + 1]
            qmax = max(qmax, z[k])
            emin = min(emin, z[k + 1])
        }
        if (z[2 * (n - 1)] < 0.0) {
            throw Exception(negZ)
        }
        d += z[2 * (n - 1)]
        if (e == 0.0) {
            for (k in 1..n) {
                z[k] = z[2 * k]
            }
            dlasrt(SortDecreasing, n, z)
            z[2 * (n - 1)] = d
            return info
        }
        let trace = d + e
        if (trace == 0.0) {
            z[2 * (n - 1)] = 0.0
            return info
        }
        for (k in 2 * n..=2 : -2) {
            z[2 * k - 1] = 0.0
            z[2 * k - 2] = z[k - 1]
            z[2 * k - 3] = 0.0
            z[2 * k - 4] = z[k - 2]
        }
        let i0 = 0
        var n0 = n - 1

        if (cbias * z[4 * i0] < z[4 * n0]) {
            let ipn4Out = 4 * (i0 + n0 + 2)
            for (i4loop in 4 * (i0 + 1)..=2 * (i0 + n0 + 1) : 4) {
                let i4 = i4loop - 1
                let ipn4 = ipn4Out - 1
                (z[i4 - 3], z[ipn4 - i4 - 4]) = (z[ipn4 - i4 - 4], z[i4 - 3])
                (z[i4 - 1], z[ipn4 - i4 - 6]) = (z[ipn4 - i4 - 6], z[i4 - 1])
            }
        }

        var pp = 0
        for (_ in 0..2) {
            d = z[4 * n0 + pp]
            for (i4loop in 4 * n0 + pp..=4 * (i0 + 1) + pp : -4) {
                let i4 = i4loop - 1
                if (z[i4 - 1] <= tol2 * d) {
                    z[i4 - 1] = copysign(0.0, -1.0)
                    d = z[i4 - 3]
                } else {
                    d = z[i4 - 3] * (d / (d + z[i4 - 1]))
                }
            }
            emin = z[4 * (i0 + 1) + pp]
            d = z[4 * i0 + pp]
            for (i4loop in 4 * (i0 + 1) + pp..=4 * n0 + pp : 4) {
                let i4 = i4loop - 1
                z[i4 - 2 * pp - 2] = d + z[i4 - 1]
                if (z[i4 - 1] <= tol2 * d) {
                    z[i4 - 1] = copysign(0.0, -1.0)
                    z[i4 - 2 * pp - 2] = d
                    z[i4 - 2 * pp] = 0.0
                    d = z[i4 + 1]
                } else if (safmin * z[i4 + 1] < z[i4 - 2 * pp - 2] && safmin * z[i4 - 2 * pp - 2] < z[i4 + 1]) {
                    let tmp = z[i4 + 1] / z[i4 - 2 * pp - 2]
                    z[i4 - 2 * pp] = z[i4 - 1] * tmp
                    d *= tmp
                } else {
                    z[i4 - 2 * pp] = z[i4 + 1] * (z[i4 - 1] / z[i4 - 2 * pp - 2])
                    d = z[i4 + 1] * (d / z[i4 - 2 * pp - 2])
                }
                emin = min(emin, z[i4 - 2 * pp])
            }
            z[4 * (n0 + 1) - pp - 3] = d

            qmax = z[4 * (i0 + 1) - pp - 3]
            for (i4loop in 4 * (i0 + 1) - pp + 2..=4 * (n0 + 1) + pp - 2 : 4) {
                let i4 = i4loop - 1
                qmax = max(qmax, z[i4])
            }
            pp = 1 - pp
        }

        var ttype: Int64 = 0
        var dmin1: Float64 = 0.0
        var dmin2: Float64 = 0.0
        var dn: Float64 = 0.0
        var dn1: Float64 = 0.0
        var dn2: Float64 = 0.0
        var g: Float64 = 0.0
        var tau: Float64 = 0.0
        var tempq: Float64
        var iter = 2
        var nFail: Int64 = 0
        var nDiv = 2 * (n0 - i0)
        var i4: Int64 = 0

        var outer = true

        for (_ in 1..=n + 1) {
            if (n0 < 0) {
                for (k in 1..n) {
                    z[k] = z[4 * k]
                }
                dlasrt(SortDecreasing, n, z)
                e = 0.0
                for (k in n - 1..=0 : -1) {
                    e += z[k]
                }
                z[2 * n] = trace
                z[2 * n + 1] = e
                z[2 * n + 2] = Float64(iter)
                z[2 * n + 3] = Float64(nDiv) / Float64(n * n)
                z[2 * n + 4] = 100.0 * Float64(nFail) / Float64(iter)
                return info
            }

            var desig: Float64 = 0.0
            var sigma: Float64 = 0.0
            if (n0 != n - 1) {
                sigma = -z[4 * (n0 + 1) - 2]
            }
            if (sigma < 0.0) {
                info = 1
                return info
            }
            var emax: Float64 = 0.0
            if (n0 > i0) {
                emin = abs(z[4 * (n0 + 1) - 6])
            } else {
                emin = 0.0
            }
            var qmin = z[4 * (n0 + 1) - 4]
            qmax = qmin
            var zSmall = false
            for (i4loop in 4 * (n0 + 1)..=8 : -4) {
                i4 = i4loop - 1
                if (z[i4 - 5] <= 0.0) {
                    zSmall = true
                    break
                }
                if (qmin >= 4.0 * emax) {
                    qmin = min(qmin, z[i4 - 3])
                    emax = max(emax, z[i4 - 5])
                }
                qmax = max(qmax, z[i4 - 7] + z[i4 - 5])
                emin = min(emin, z[i4 - 5])
            }
            if (!zSmall) {
                i4 = 3
            }
            var i0 = (i4 + 1) / 4 - 1
            pp = 0
            if (n0 - i0 > 1) {
                var dee = z[4 * i0]
                var deemin = dee
                var kmin = i0
                for (i4loop in 4 * (i0 + 1) + 1..=4 * (n0 + 1) - 3 : 4) {
                    i4 = i4loop - 1
                    dee = z[i4] * (dee / (dee + z[i4 - 2]))
                    if (dee <= deemin) {
                        deemin = dee
                        kmin = (i4 + 4) / 4 - 1
                    }
                }
                if ((kmin - i0) * 2 < n0 - kmin && deemin <= 0.5 * z[4 * n0]) {
                    let ipn4Out = 4 * (i0 + n0 + 2)
                    pp = 2
                    for (i4loop in 4 * (i0 + 1)..=2 * (i0 + n0 + 1) : 4) {
                        i4 = i4loop - 1
                        let ipn4 = ipn4Out - 1
                        (z[i4 - 3], z[ipn4 - i4 - 4]) = (z[ipn4 - i4 - 4], z[i4 - 3])
                        (z[i4 - 2], z[ipn4 - i4 - 3]) = (z[ipn4 - i4 - 3], z[i4 - 2])
                        (z[i4 - 1], z[ipn4 - i4 - 6]) = (z[ipn4 - i4 - 6], z[i4 - 1])
                        (z[i4], z[ipn4 - i4 - 5]) = (z[ipn4 - i4 - 5], z[i4])
                    }
                }
            }
            var dmin: Float64 = Float64(-max(0, Int64(qmin - 2.0 * sqrt(qmin) * sqrt(emax))))

            var nbig = 100 * (n0 - i0 + 1)
            for (_ in 0..nbig) {
                if (i0 > n0) {
                    outer = false
                    break
                }

                (i0, n0, pp, dmin, sigma, desig, qmax, nFail, iter, nDiv, ttype, dmin1, dmin2, dn, dn1, dn2, g, tau) = dlasq3(
                    i0, n0, z, pp, dmin, sigma, desig, qmax, nFail, iter, nDiv, ttype, dmin1, dmin2, dn, dn1, dn2, g,
                    tau)

                pp = 1 - pp
                if (pp == 0 && n0 - i0 >= 3) {
                    if (z[4 * (n0 + 1) - 1] <= tol2 * qmax || z[4 * (n0 + 1) - 2] <= tol2 * sigma) {
                        var splt = i0 - 1
                        qmax = z[4 * i0]
                        emin = z[4 * (i0 + 1) - 2]
                        var oldemn = z[4 * (i0 + 1) - 1]
                        for (i4loop in 4 * (i0 + 1)..=4 * (n0 - 2) : 4) {
                            i4 = i4loop - 1
                            if (z[i4] <= tol2 * z[i4 - 3] || z[i4 - 1] <= tol2 * sigma) {
                                z[i4 - 1] = -sigma
                                splt = i4 / 4
                                qmax = 0.0
                                emin = z[i4 + 3]
                                oldemn = z[i4 + 4]
                            } else {
                                qmax = max(qmax, z[i4 + 1])
                                emin = min(emin, z[i4 - 1])
                                oldemn = min(oldemn, z[i4])
                            }
                        }
                        z[4 * (n0 + 1) - 2] = emin
                        z[4 * (n0 + 1) - 1] = oldemn
                        i0 = splt + 1
                    }
                }
            }

            if (!outer) {
                continue
            }

            info = 2
            i1 = i0
            while (true) {
                tempq = z[4 * i0]
                z[4 * i0] += sigma
                for (k in i0 + 1..=n0) {
                    let tempe = z[4 * (k + 1) - 6]
                    z[4 * (k + 1) - 6] *= tempq / z[4 * (k + 1) - 8]
                    tempq = z[4 * k]
                    z[4 * k] += sigma + tempe - z[4 * (k + 1) - 6]
                }
                if (i1 <= 0) {
                    break
                }
                n1 = i1 - 1
                while (i1 >= 1 && z[4 * (i1 + 1) - 6] >= 0.0) {
                    i1 -= 1
                }
                sigma = -z[4 * (n1 + 1) - 2]
            }
            for (k in 0..n) {
                z[2 * k] = z[4 * k]
                if (k < n0) {
                    z[2 * (k + 1) - 1] = z[4 * (k + 1) - 1]
                } else {
                    z[2 * (k + 1)] = 0.0
                }
            }
            return info
        }
        info = 3
        return info
    }
}
