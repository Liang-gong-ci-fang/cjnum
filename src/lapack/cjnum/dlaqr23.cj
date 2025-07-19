package cjnum.lapack.cjnum

import cjnum.blas.blas64.*
import cjnum.blas.*
import cjnum.lapack.*
import std.math.*

extend LFloat64Implementation {
    public func dlaqr23(wantt: Bool, wantz: Bool, n: Int64, ktop: Int64, kbot: Int64, nw: Int64, h: Array<Float64>,
        ldh: Int64, iloz: Int64, ihiz: Int64, z: Array<Float64>, ldz: Int64, sr: Array<Float64>, si: Array<Float64>,
        v: Array<Float64>, ldv: Int64, nh: Int64, t: Array<Float64>, ldt: Int64, nv: Int64, wv: Array<Float64>,
        ldwv: Int64, work: Array<Float64>, lwork: Int64, recur: Int64): (Int64, Int64) {
        var ns: Int64 = 0
        var nd: Int64 = 0

        if (n < 0) {
            throw Exception(nLT0)
        } else if (ktop < 0 || max(0, n - 1) < ktop) {
            throw Exception(badKtop)
        } else if (kbot < min(ktop, n - 1) || n <= kbot) {
            throw Exception(badKbot)
        } else if (nw < 0 || kbot - ktop + 1 + 1 < nw) {
            throw Exception(badNw)
        } else if (ldh < max(1, n)) {
            throw Exception(badLdH)
        } else if (wantz && (iloz < 0 || ktop < iloz)) {
            throw Exception(badIloz)
        } else if (wantz && (ihiz < kbot || n <= ihiz)) {
            throw Exception(badIhiz)
        } else if (ldz < 1 || (wantz && ldz < n)) {
            throw Exception(badLdZ)
        } else if (ldv < max(1, nw)) {
            throw Exception(badLdV)
        } else if (nh < nw) {
            throw Exception(badNh)
        } else if (ldt < max(1, nh)) {
            throw Exception(badLdT)
        } else if (nv < 0) {
            throw Exception(nvLT0)
        } else if (ldwv < max(1, nw)) {
            throw Exception(badLdWV)
        } else if (lwork < max(1, 2 * nw) && lwork != -1) {
            throw Exception(badLWork)
        } else if (work.size < max(1, lwork)) {
            throw Exception(shortWork)
        } else if (recur < 0) {
            throw Exception(recurLT0)
        }

        if (nw == 0) {
            work[0] = 1.0
            return (0, 0)
        }

        let jw = nw
        var lwkopt = max(1, 2 * nw)
        if (jw > 2) {
            dgehrd(jw, 0, jw - 2, t, ldt, work, work, -1)
            let lwk1 = Int64(work[0])

            dormhr(Right, NoTrans, jw, jw, 0, jw - 2, t, ldt, work, v, ldv, work, -1)
            let lwk2 = Int64(work[0])
            if (recur > 0) {
                dlaqr04(true, true, jw, 0, jw - 1, t, ldt, sr, si, 0, jw - 1, v, ldv, work, -1, recur - 1)
                let lwk3 = Int64(work[0])
                lwkopt = max(jw + max(lwk1, lwk2), lwk3)
            } else {
                lwkopt = jw + max(lwk1, lwk2)
            }
        }

        if (lwork == -1) {
            work[0] = Float64(lwkopt)
            return (0, 0)
        }

        if (h.size < (n - 1) * ldh + n) {
            throw Exception(shortH)
        } else if (v.size < (nw - 1) * ldv + nw) {
            throw Exception(shortV)
        } else if (t.size < (nw - 1) * ldt + nh) {
            throw Exception(shortT)
        } else if (wv.size < (nv - 1) * ldwv + nw) {
            throw Exception(shortWV)
        } else if (wantz && z.size < (n - 1) * ldz + n) {
            throw Exception(shortZ)
        } else if (sr.size != kbot + 1) {
            throw Exception(badLenSr)
        } else if (si.size != kbot + 1) {
            throw Exception(badLenSi)
        } else if (ktop > 0 && h[ktop * ldh + ktop - 1] != 0.0) {
            throw Exception(notIsolated)
        } else if (kbot + 1 < n && h[(kbot + 1) * ldh + kbot] != 0.0) {
            throw Exception(notIsolated)
        }

        let ulp = dlamchP
        let smlnum = Float64(n) / ulp * dlamchS

        var s = 0.0
        let kwtop = kbot - jw + 1
        if (kwtop != ktop) {
            s = h[kwtop * ldh + kwtop - 1]
        }
        if (kwtop == kbot) {
            sr[kwtop] = h[kwtop * ldh + kwtop]
            si[kwtop] = 0.0
            ns = 1
            nd = 0
            if (abs(s) <= max(smlnum, ulp * abs(h[kwtop * ldh + kwtop]))) {
                ns = 0
                nd = 1
                if (kwtop > ktop) {
                    h[kwtop * ldh + kwtop - 1] = 0.0
                }
            }
            work[0] = 1.0
            return (ns, nd)
        }

        dlacpy(Upper, jw, jw, h[kwtop * ldh + kwtop..], ldh, t, ldt)
        let nFloat64 = nFloat64Implementation()
        nFloat64.dcopy(jw - 1, h[(kwtop + 1) * ldh + kwtop..], ldh + 1, t[ldt..], ldt + 1)
        dlaset(All, jw, jw, 0.0, 1.0, v, ldv)
        let nmin = ilaenv(12, "DLAQR3", "SV", jw, 0, jw - 1, lwork)
        var infqr = 0
        if (recur > 0 && jw > nmin) {
            infqr = dlaqr04(true, true, jw, 0, jw - 1, t, ldt, sr[kwtop..], si[kwtop..], 0, jw - 1, v, ldv, work, lwork,
                recur - 1)
        } else {
            infqr = dlahqr(true, true, jw, 0, jw - 1, t, ldt, sr[kwtop..], si[kwtop..], 0, jw - 1, v, ldv)
        }

        for (j in 0..jw - 3) {
            t[(j + 2) * ldt + j] = 0.0
            t[(j + 3) * ldt + j] = 0.0
        }
        if (jw >= 3) {
            t[(jw - 1) * ldt + jw - 3] = 0.0
        }

        ns = jw
        var ilst = infqr

        while (ilst < ns) {
            var bulge = false
            if (ns >= 2) {
                bulge = t[(ns - 1) * ldt + ns - 2] != 0.0
            }
            if (!bulge) {
                var abst = abs(t[(ns - 1) * ldt + ns - 1])
                if (abst == 0.0) {
                    abst = abs(s)
                }
                if (abs(s * v[ns - 1]) <= max(smlnum, ulp * abst)) {
                    ns--
                } else {
                    (_, ilst, _) = dtrexc(UpdateSchur, jw, t, ldt, v, ldv, ns - 1, ilst, work)
                    ilst++
                }
                continue
            }
            var abst = abs(t[(ns - 1) * ldt + ns - 1]) + sqrt(abs(t[(ns - 1) * ldt + ns - 2])) * sqrt(
                abs(t[(ns - 2) * ldt + ns - 1]))
            if (abst == 0.0) {
                abst = abs(s)
            }
            if (max(abs(s * v[ns - 1]), abs(s * v[ns - 2])) <= max(smlnum, ulp * abst)) {
                ns -= 2
            } else {
                (_, ilst, _) = dtrexc(UpdateSchur, jw, t, ldt, v, ldv, ns - 1, ilst, work)
                ilst += 2
            }
        }

        if (ns == 0) {
            s = 0.0
        }
        if (ns < jw) {
            var sorted = false
            var i = ns
            while (!sorted) {
                sorted = true
                let kend = i - 1
                i = infqr
                var k = 0
                if (i == ns - 1 || t[(i + 1) * ldt + i] == 0.0) {
                    k = i + 1
                } else {
                    k = i + 2
                }
                while (k <= kend) {
                    var evi = 0.0
                    if (k == i + 1) {
                        evi = abs(t[i * ldt + i])
                    } else {
                        evi = abs(t[i * ldt + i]) + sqrt(abs(t[(i + 1) * ldt + i])) * sqrt(abs(t[i * ldt + i + 1]))
                    }

                    var evk = 0.0
                    if (k == kend || t[(k + 1) * ldt + k] == 0.0) {
                        evk = abs(t[k * ldt + k])
                    } else {
                        evk = abs(t[k * ldt + k]) + sqrt(abs(t[(k + 1) * ldt + k])) * sqrt(abs(t[k * ldt + k + 1]))
                    }

                    if (evi >= evk) {
                        i = k
                    } else {
                        sorted = false
                        let (_, ilst, ok) = dtrexc(UpdateSchur, jw, t, ldt, v, ldv, i, k, work)
                        if (ok) {
                            i = ilst
                        } else {
                            i = k
                        }
                    }
                    if (i == kend || t[(i + 1) * ldt + i] == 0.0) {
                        k = i + 1
                    } else {
                        k = i + 2
                    }
                }
            }
        }

        var i = jw - 1
        while (i >= infqr) {
            if (i == infqr || t[i * ldt + i - 1] == 0.0) {
                sr[kwtop + i] = t[i * ldt + i]
                si[kwtop + i] = 0.0
                i--
                continue
            }
            let aa = t[(i - 1) * ldt + i - 1]
            let bb = t[(i - 1) * ldt + i]
            let cc = t[i * ldt + i - 1]
            let dd = t[i * ldt + i]
            (_, _, _, _, sr[kwtop + i - 1], si[kwtop + i - 1], sr[kwtop + i], si[kwtop + i], _, _) = dlanv2(aa, bb, cc,
                dd)
            i -= 2
        }

        if (ns < jw || s == 0.0) {
            if (ns > 1 && s != 0.0) {
                nFloat64.dcopy(ns, v[0..ns], 1, work[0..ns], 1)
                let (_, tau) = dlarfg(ns, work[0], work[1..ns], 1)
                work[0] = 1.0
                dlaset(Lower, jw - 2, jw - 2, 0.0, 0.0, t[2 * ldt..], ldt)
                dlarf(Left, ns, jw, work[0..ns], 1, tau, t, ldt, work[jw..])
                dlarf(Right, ns, ns, work[0..ns], 1, tau, t, ldt, work[jw..])
                dlarf(Right, jw, ns, work[0..ns], 1, tau, v, ldv, work[jw..])
                dgehrd(jw, 0, ns - 1, t, ldt, work[0..jw - 1], work[jw..], lwork - jw)
            }

            if (kwtop > 0) {
                h[kwtop * ldh + kwtop - 1] = s * v[0]
            }
            dlacpy(Upper, jw, jw, t, ldt, h[kwtop * ldh + kwtop..], ldh)
            nFloat64.dcopy(jw - 1, t[ldt..], ldt + 1, h[(kwtop + 1) * ldh + kwtop..], ldh + 1)

            if (ns > 1 && s != 0.0) {
                dormhr(Right, NoTrans, jw, ns, 0, ns - 1, t, ldt, work[0..ns - 1], v, ldv, work[jw..], lwork - jw)
            }

            var ltop = 0
            if (!wantt) {
                ltop = ktop
            }
            var krow = ltop
            while (krow < kwtop) {
                let kln = min(nv, kwtop - krow)
                nFloat64.dgemm(NoTrans, NoTrans, kln, jw, jw, 1.0, h[krow * ldh + kwtop..], ldh, v, ldv, 0.0, wv, ldwv)
                dlacpy(All, kln, jw, wv, ldwv, h[krow * ldh + kwtop..], ldh)
                krow += nv
            }

            if (wantt) {
                var kcol = kbot + 1
                while (kcol < n) {
                    let kln = min(nh, n - kcol)
                    nFloat64.dgemm(Trans, NoTrans, jw, kln, jw, 1.0, v, ldv, h[kwtop * ldh + kcol..], ldh, 0.0, t, ldt)
                    dlacpy(All, jw, kln, t, ldt, h[kwtop * ldh + kcol..], ldh)
                    kcol += nh
                }
            }

            if (wantz) {
                var krow = iloz
                while (krow <= ihiz) {
                    let kln = min(nv, ihiz - krow + 1)
                    nFloat64.dgemm(NoTrans, NoTrans, kln, jw, jw, 1.0, z[krow * ldz + kwtop..], ldz, v, ldv, 0.0, wv,
                        ldwv)
                    dlacpy(All, kln, jw, wv, ldwv, z[krow * ldz + kwtop..], ldz)
                    krow += nv
                }
            }
        }

        nd = jw - ns

        ns -= infqr
        work[0] = Float64(lwkopt)
        return (ns, nd)
    }
}
