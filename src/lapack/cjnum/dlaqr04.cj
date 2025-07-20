package cjnum.lapack.cjnum

import cjnum.blas.*
import std.math.*

extend LFloat64Implementation {
    public func dlaqr04(wantt: Bool, wantz: Bool, n: Int64, ilo: Int64, ihi: Int64, h: Array<Float64>, ldh: Int64,
        wr: Array<Float64>, wi: Array<Float64>, iloz: Int64, ihiz: Int64, z: Array<Float64>, ldz: Int64,
        work: Array<Float64>, lwork: Int64, recur: Int64): Int64 {
        var unconverged: Int64 = 0
        const ntiny = 15
        const kexnw = 5
        const kexsh = 6
        const wilk1 = 0.75
        const wilk2 = -0.4375

        if (n < 0) {
            throw Exception(nLT0)
        } else if ((ilo < 0) || (max(0, n - 1) < ilo)) {
            throw Exception(badIlo)
        } else if ((ihi < min(ilo, n - 1)) || (n <= ihi)) {
            throw Exception(badIhi)
        } else if (ldh < max(1, n)) {
            throw Exception(badLdH)
        } else if (wantz && ((iloz < 0) || (ilo < iloz))) {
            throw Exception(badIloz)
        } else if (wantz && ((ihiz < ihi) || (n <= ihiz))) {
            throw Exception(badIhiz)
        } else if ((ldz < 1) || (wantz && ldz < n)) {
            throw Exception(badLdZ)
        } else if ((lwork < 1) && (lwork != -1)) {
            throw Exception(badLWork)
        } else if (work.size < max(1, lwork)) {
            throw Exception(shortWork)
        } else if (recur < 0) {
            throw Exception(recurLT0)
        }

        if (n == 0) {
            work[0] = 1.0
            return 0
        }

        if (lwork != -1) {
            if (h.size < (n - 1) * ldh + n) {
                throw Exception(shortH)
            }
            if (wr.size != ihi + 1) {
                throw Exception(badLenWr)
            }
            if (wi.size != ihi + 1) {
                throw Exception(badLenWi)
            }
            if (wantz && z.size < (n - 1) * ldz + n) {
                throw Exception(shortZ)
            }
            if (ilo > 0 && h[ilo * ldh + ilo - 1] != 0.0) {
                throw Exception(notIsolated)
            }
            if (ihi + 1 < n && h[(ihi + 1) * ldh + ihi] != 0.0) {
                throw Exception(notIsolated)
            }
        }

        if (n <= ntiny) {
            if (lwork == -1) {
                work[0] = 1.0
                return 0
            }
            return dlahqr(wantt, wantz, n, ilo, ihi, h, ldh, wr, wi, iloz, ihiz, z, ldz)
        }

        var jbcmpz: String
        if (wantt) {
            jbcmpz = "S"
        } else {
            jbcmpz = "E"
        }
        if (wantz) {
            jbcmpz += "V"
        } else {
            jbcmpz += "N"
        }

        var fname: String
        if (recur > 0) {
            fname = "DLAQR0"
        } else {
            fname = "DLAQR4"
        }

        var nwr = ilaenv(13, fname, jbcmpz, n, ilo, ihi, lwork)
        nwr = max(2, nwr)
        nwr = min(ihi - ilo + 1, min((n - 1) / 3, nwr))

        var nsr = ilaenv(15, fname, jbcmpz, n, ilo, ihi, lwork)
        nsr = min(nsr, min((n - 3) / 6, ihi - ilo))
        nsr = max(2, nsr & !1)

        dlaqr23(wantt, wantz, n, ilo, ihi, nwr + 1, h, ldh, iloz, ihiz, z, ldz, wr, wi, h, ldh, n, h, ldh, n, h, ldh,
            work, -1, recur)
        let lwkopt = max(3 * nsr / 2, Int64(work[0]))
        if (lwork == -1) {
            work[0] = Float64(lwkopt)
            return 0
        }

        var nmin = ilaenv(12, fname, jbcmpz, n, ilo, ihi, lwork)
        nmin = max(ntiny, nmin)

        var nibble = ilaenv(14, fname, jbcmpz, n, ilo, ihi, lwork)
        nibble = max(0, nibble)

        var kacc22 = ilaenv(16, fname, jbcmpz, n, ilo, ihi, lwork)
        kacc22 = max(0, min(kacc22, 2))

        let nwmax = min((n - 1) / 3, lwork / 2)
        var nw = nwmax

        let nsmax = min((n - 3) / 6, 2 * lwork / 3) & !1

        var ndfl = 1
        var ndec = 0

        let itmax = max(30, 2 * kexsh) * max(10, (ihi - ilo + 1))
        var it = 0
        var kbot = ihi
        while (kbot >= ilo) {
            if (it == itmax) {
                unconverged = kbot + 1
                break
            }
            it++

            var ktop = ilo
            for (k in kbot..=ilo + 1 : -1) {
                if (h[k * ldh + k - 1] == 0.0) {
                    ktop = k
                    break
                }
            }

            let nh = kbot - ktop + 1
            let nwupbd = min(nh, nwmax)
            if (ndfl < kexnw) {
                nw = min(nwupbd, nwr)
            } else {
                nw = min(nwupbd, 2 * nw)
            }
            if (nw < nwmax) {
                if (nw >= nh - 1) {
                    nw = nh
                } else {
                    let kwtop = kbot - nw + 1
                    if (abs(h[kwtop * ldh + kwtop - 1]) > abs(h[(kwtop - 1) * ldh + kwtop - 2])) {
                        nw++
                    }
                }
            }
            if (ndfl < kexnw) {
                ndec = -1
            } else if (ndec >= 0 || nw >= nwupbd) {
                ndec++
                if (nw - ndec < 2) {
                    ndec = 0
                }
                nw -= ndec
            }

            let kv = n - nw
            var kt = nw
            var kwv = nw + 1
            var nhv = n - kwv - kt

            let (ls, ld) = dlaqr23(wantt, wantz, n, ktop, kbot, nw, h, ldh, iloz, ihiz, z, ldz, wr[..kbot + 1],
                wi[..kbot + 1], h[kv * ldh..], ldh, nhv, h[kv * ldh + kt..], ldh, nhv, h[kwv * ldh..], ldh, work, lwork,
                recur)

            kbot -= ld
            var ks = kbot - ls + 1

            if (ld > 0 && (100 * ld > nw * nibble || kbot - ktop + 1 <= min(nmin, nwmax))) {
                ndfl = 1
                continue
            }

            var ns = min(min(nsmax, nsr), max(2, kbot - ktop)) & !1

            if (ndfl % kexsh == 0) {
                ks = kbot - ns + 1
                for (i in kbot..max(ks, ktop + 1) : -2) {
                    let ss = abs(h[i * ldh + i - 1]) + abs(h[(i - 1) * ldh + i - 2])
                    let aa = wilk1 * ss + h[i * ldh + i]
                    (_, _, _, _, wr[i - 1], wi[i - 1], wr[i], wi[i], _, _) = dlanv2(aa, ss, wilk2 * ss, aa)
                }
                if (ks == ktop) {
                    wr[ks + 1] = h[(ks + 1) * ldh + ks + 1]
                    wi[ks + 1] = 0.0
                    wr[ks] = wr[ks + 1]
                    wi[ks] = wi[ks + 1]
                }
            } else {
                if (kbot - ks + 1 <= ns / 2) {
                    ks = kbot - ns + 1
                    kt = n - ns
                    dlacpy(All, ns, ns, h[ks * ldh + ks..], ldh, h[kt * ldh..], ldh)
                    if (ns > nmin && recur > 0) {
                        ks += dlaqr04(false, false, ns, 1, ns - 1, h[kt * ldh..], ldh, wr[ks..ks + ns], wi[ks..ks + ns],
                            0, 0, [], 0, work, lwork, recur - 1)
                    } else {
                        ks += dlahqr(false, false, ns, 0, ns - 1, h[kt * ldh..], ldh, wr[ks..ks + ns], wi[ks..ks + ns],
                            0, 0, [], 1)
                    }
                    if (ks >= kbot) {
                        let aa = h[(kbot - 1) * ldh + kbot - 1]
                        let bb = h[(kbot - 1) * ldh + kbot]
                        let cc = h[kbot * ldh + kbot - 1]
                        let dd = h[kbot * ldh + kbot]
                        (_, _, _, _, wr[kbot - 1], wi[kbot - 1], wr[kbot], wi[kbot], _, _) = dlanv2(aa, bb, cc, dd)
                        ks = kbot - 1
                    }
                }

                if (kbot - ks + 1 > ns) {
                    var sorted = false
                    for (k in kbot..ks : -1) {
                        if (sorted) {
                            break
                        }
                        sorted = true
                        for (i in ks..k) {
                            if (abs(wr[i]) + abs(wi[i]) >= abs(wr[i + 1]) + abs(wi[i + 1])) {
                                continue
                            }
                            sorted = false
                            (wr[i], wr[i + 1]) = (wr[i + 1], wr[i])
                            (wi[i], wi[i + 1]) = (wi[i + 1], wi[i])
                        }
                    }
                }

                for (i in kbot..ks + 1 : -2) {
                    if (wi[i] == -wi[i - 1]) {
                        continue
                    }
                    (wr[i], wr[i - 1], wr[i - 2]) = (wr[i - 1], wr[i - 2], wr[i])
                    (wi[i], wi[i - 1], wi[i - 2]) = (wi[i - 1], wi[i - 2], wi[i])
                }
            }

            if (kbot - ks + 1 == 2 && wi[kbot] == 0.0) {
                if (abs(wr[kbot] - h[kbot * ldh + kbot]) < abs(wr[kbot - 1] - h[kbot * ldh + kbot])) {
                    wr[kbot - 1] = wr[kbot]
                } else {
                    wr[kbot] = wr[kbot - 1]
                }
            }

            ns = min(ns, kbot - ks + 1) & !1
            ks = kbot - ns + 1

            let kdu = 2 * ns
            let ku = n - kdu
            let kwh = kdu
            kwv = kdu + 3
            nhv = n - kwv - kdu

            dlaqr5(wantt, wantz, kacc22, n, ktop, kbot, ns, wr[ks..ks + ns], wi[ks..ks + ns], h, ldh, iloz, ihiz, z, ldz,
                work, 3, h[ku * ldh..], ldh, nhv, h[kwv * ldh..], ldh, nhv, h[ku * ldh + kwh..], ldh)

            if (ld > 0) {
                ndfl = 1
            } else {
                ndfl++
            }
        }

        work[0] = Float64(lwkopt)
        return unconverged
    }
}
