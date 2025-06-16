package cjnum.lapack.cjnum

import cjnum.blas.blas64.*
import std.math.*

extend LFloat64Implementation {
    public func dlahqr(wantt: Bool, wantz: Bool, n: Int64, ilo: Int64, ihi: Int64, h: Array<Float64>, ldh: Int64,
        wr: Array<Float64>, wi: Array<Float64>, iloz: Int64, ihiz: Int64, z: Array<Float64>, ldz: Int64): Int64 {
        if (n < 0) {
            throw Exception(nLT0)
        } else if ((ilo < 0) || (max(0, ihi) < ilo)) {
            throw Exception(badIlo)
        } else if (ihi >= n) {
            throw Exception(badIhi)
        } else if (ldh < max(1, n)) {
            throw Exception(badLdH)
        } else if (wantz && ((iloz < 0) || (ilo < iloz))) {
            throw Exception(badIloz)
        } else if (wantz && ((ihiz < ihi) || (n <= ihiz))) {
            throw Exception(badIhiz)
        } else if ((ldz < 1) || (wantz && ldz < n)) {
            throw Exception(badLdZ)
        }

        if (n == 0) {
            return 0
        }

        if (h.size < (n - 1) * ldh + n) {
            throw Exception(shortH)
        } else if (wr.size != ihi + 1) {
            throw Exception(shortWr)
        } else if (wi.size != ihi + 1) {
            throw Exception(shortWi)
        } else if (wantz && z.size < (n - 1) * ldz + n) {
            throw Exception(shortZ)
        } else if ((ilo > 0) && (h[ilo * ldh + ilo - 1] != 0.0)) {
            throw Exception(notIsolated)
        }

        if (ilo == ihi) {
            wr[ilo] = h[ilo * ldh + ilo]
            wi[ilo] = 0.0
            return 0
        }

        for (j in ilo..ihi - 2) {
            h[(j + 2) * ldh + j] = 0.0
            h[(j + 3) * ldh + j] = 0.0
        }
        if (ilo <= ihi - 2) {
            h[ihi * ldh + ihi - 2] = 0.0
        }

        let nh = ihi - ilo + 1
        let nz = ihiz - iloz + 1

        let ulp = dlamchP
        let smlnum = Float64(nh) / ulp * dlamchS

        var i1: Int64 = 0
        var i2: Int64 = 0
        if (wantt) {
            i1 = 0
            i2 = n - 1
        }

        let itmax = 30 * max(10, nh)
        var kdefl = 0

        let nFloat64 = nFloat64Implementation()
        var i = ihi
        while (i >= ilo) {
            var l = ilo

            var converged = false
            for (_ in 0..=itmax) {
                var k: Int64 = i
                while (k > l) {
                    if (abs(h[k * ldh + k - 1]) <= smlnum) {
                        break
                    }
                    var tst = abs(h[(k - 1) * ldh + k - 1]) + abs(h[k * ldh + k])
                    if (tst == 0.0) {
                        if (k - 2 >= ilo) {
                            tst += abs(h[(k - 1) * ldh + k - 2])
                        }
                        if (k + 1 <= ihi) {
                            tst += abs(h[(k + 1) * ldh + k])
                        }
                    }
                    if (abs(h[k * ldh + k - 1]) <= ulp * tst) {
                        let ab = max(abs(h[k * ldh + k - 1]), abs(h[(k - 1) * ldh + k]))
                        let ba = min(abs(h[k * ldh + k - 1]), abs(h[(k - 1) * ldh + k]))
                        let aa = max(abs(h[k * ldh + k]), abs(h[(k - 1) * ldh + k - 1] - h[k * ldh + k]))
                        let bb = min(abs(h[k * ldh + k]), abs(h[(k - 1) * ldh + k - 1] - h[k * ldh + k]))
                        let s = aa + ab
                        if (ab / s * ba <= max(smlnum, aa / s * bb * ulp)) {
                            break
                        }
                    }
                    k--
                }
                l = k
                if (l > ilo) {
                    h[l * ldh + l - 1] = 0.0
                }
                if (l >= i - 1) {
                    converged = true
                    break
                }
                kdefl++

                if (!wantt) {
                    i1 = l
                    i2 = i
                }

                const dat1 = 0.75
                const dat2 = -0.4375
                const kexsh = 10
                var h11: Float64 = 0.0
                var h21: Float64 = 0.0
                var h12: Float64 = 0.0
                var h22: Float64 = 0.0
                if (kdefl % (2 * kexsh) == 0) {
                    let s = abs(h[i * ldh + i - 1]) + abs(h[(i - 1) * ldh + i - 2])
                    h11 = dat1 * s + h[i * ldh + i]
                    h12 = dat2 * s
                    h21 = s
                    h22 = h11
                } else if (kdefl % kexsh == 0) {
                    let s = abs(h[(l + 1) * ldh + l]) + abs(h[(l + 2) * ldh + l + 1])
                    h11 = dat1 * s + h[l * ldh + l]
                    h12 = dat2 * s
                    h21 = s
                    h22 = h11
                } else {
                    h11 = h[(i - 1) * ldh + i - 1]
                    h21 = h[i * ldh + i - 1]
                    h12 = h[(i - 1) * ldh + i]
                    h22 = h[i * ldh + i]
                }
                let s = abs(h11) + abs(h12) + abs(h21) + abs(h22)
                var rt1r: Float64 = 0.0
                var rt1i: Float64 = 0.0
                var rt2r: Float64 = 0.0
                var rt2i: Float64 = 0.0
                if (s != 0.0) {
                    h11 /= s
                    h21 /= s
                    h12 /= s
                    h22 /= s
                    let tr = (h11 + h22) / 2.0
                    let det = (h11 - tr) * (h22 - tr) - h12 * h21
                    let rtdisc = sqrt(abs(det))
                    if (det >= 0.0) {
                        rt1r = tr * s
                        rt2r = rt1r
                        rt1i = rtdisc * s
                        rt2i = -rt1i
                    } else {
                        rt1r = tr + rtdisc
                        rt2r = tr - rtdisc
                        if (abs(rt1r - h22) <= abs(rt2r - h22)) {
                            rt1r *= s
                            rt2r = rt1r
                        } else {
                            rt2r *= s
                            rt1r = rt2r
                        }
                        rt1i = 0.0
                        rt2i = 0.0
                    }
                }

                var m: Int64 = i - 2
                let v = Array<Float64>(3, item: 0.0)
                while (m >= l) {
                    var h21s = h[(m + 1) * ldh + m]
                    var s = abs(h[m * ldh + m] - rt2r) + abs(rt2i) + abs(h21s)
                    h21s /= s
                    v[0] = h21s * h[m * ldh + m + 1] + (h[m * ldh + m] - rt1r) * ((h[m * ldh + m] - rt2r) / s) - (rt2i /
                        s * rt1i)
                    v[1] = h21s * (h[m * ldh + m] + h[(m + 1) * ldh + m + 1] - rt1r - rt2r)
                    v[2] = h21s * h[(m + 2) * ldh + m + 1]
                    s = abs(v[0]) + abs(v[1]) + abs(v[2])
                    v[0] /= s
                    v[1] /= s
                    v[2] /= s
                    if (m == l) {
                        break
                    }
                    let dsum = abs(h[(m - 1) * ldh + m - 1]) + abs(h[m * ldh + m]) + abs(h[(m + 1) * ldh + m + 1])
                    if (abs(h[m * ldh + m - 1]) * (abs(v[1]) + abs(v[2])) <= ulp * abs(v[0]) * dsum) {
                        break
                    }
                    m--
                }

                for (k in m..i) {
                    let nr = min(3, i - k + 1)
                    if (k > m) {
                        nFloat64.dcopy(nr, h[k * ldh + k - 1..], ldh, v[..], 1)
                    }
                    var t0: Float64 = 0.0
                    (v[0], t0) = dlarfg(nr, v[0], v[1..], 1)
                    if (k > m) {
                        h[k * ldh + k - 1] = v[0]
                        h[(k + 1) * ldh + k - 1] = 0.0
                        if (k < i - 1) {
                            h[(k + 2) * ldh + k - 1] = 0.0
                        }
                    } else if (m > l) {
                        h[k * ldh + k - 1] *= 1.0 - t0
                    }

                    let t1 = t0 * v[1]

                    if (nr == 3) {
                        let t2 = t0 * v[2]

                        for (j in k..=i2) {
                            let sum = h[k * ldh + j] + v[1] * h[(k + 1) * ldh + j] + v[2] * h[(k + 2) * ldh + j]
                            h[k * ldh + j] -= sum * t0
                            h[(k + 1) * ldh + j] -= sum * t1
                            h[(k + 2) * ldh + j] -= sum * t2
                        }

                        for (j in i1..=min(k + 3, i)) {
                            let sum = h[j * ldh + k] + v[1] * h[j * ldh + k + 1] + v[2] * h[j * ldh + k + 2]
                            h[j * ldh + k] -= sum * t0
                            h[j * ldh + k + 1] -= sum * t1
                            h[j * ldh + k + 2] -= sum * t2
                        }

                        if (wantz) {
                            for (j in iloz..=ihiz) {
                                let sum = z[j * ldz + k] + v[1] * z[j * ldz + k + 1] + v[2] * z[j * ldz + k + 2]
                                z[j * ldz + k] -= sum * t0
                                z[j * ldz + k + 1] -= sum * t1
                                z[j * ldz + k + 2] -= sum * t2
                            }
                        }
                    } else if (nr == 2) {
                        for (j in k..=i2) {
                            let sum = h[k * ldh + j] + v[1] * h[(k + 1) * ldh + j]
                            h[k * ldh + j] -= sum * t0
                            h[(k + 1) * ldh + j] -= sum * t1
                        }

                        for (j in i1..=i) {
                            let sum = h[j * ldh + k] + v[1] * h[j * ldh + k + 1]
                            h[j * ldh + k] -= sum * t0
                            h[j * ldh + k + 1] -= sum * t1
                        }

                        if (wantz) {
                            for (j in iloz..=ihiz) {
                                let sum = z[j * ldz + k] + v[1] * z[j * ldz + k + 1]
                                z[j * ldz + k] -= sum * t0
                                z[j * ldz + k + 1] -= sum * t1
                            }
                        }
                    }
                }
            }

            if (!converged) {
                return i + 1
            }

            if (l == i) {
                wr[i] = h[i * ldh + i]
                wi[i] = 0.0
            } else if (l == i - 1) {
                var cs: Float64
                var sn: Float64
                var (a, b) = (h[(i - 1) * ldh + i - 1], h[(i - 1) * ldh + i])
                var (c, d) = (h[i * ldh + i - 1], h[i * ldh + i])
                (a, b, c, d, wr[i - 1], wi[i - 1], wr[i], wi[i], cs, sn) = dlanv2(a, b, c, d)
                (h[(i - 1) * ldh + i - 1], h[(i - 1) * ldh + i]) = (a, b)
                (h[i * ldh + i - 1], h[i * ldh + i]) = (c, d)

                if (wantt) {
                    if (i2 > i) {
                        nFloat64.drot(i2 - i, h[(i - 1) * ldh + i + 1..], 1, h[i * ldh + i + 1..], 1, cs, sn)
                    }
                    nFloat64.drot(i - i1 - 1, h[i1 * ldh + i - 1..], ldh, h[i1 * ldh + i..], ldh, cs, sn)
                }

                if (wantz) {
                    nFloat64.drot(nz, z[iloz * ldz + i - 1..], ldz, z[iloz * ldz + i..], ldz, cs, sn)
                }
            }

            kdefl = 0

            i = l - 1
        }
        return 0
    }
}
