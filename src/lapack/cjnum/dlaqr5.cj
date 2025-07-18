package cjnum.lapack.cjnum

import cjnum.blas.blas64.*
import cjnum.blas.*
import cjnum.lapack.*
import std.math.*

extend LFloat64Implementation {
    public func dlaqr5(wantt: Bool, wantz: Bool, kacc22: Int64, n: Int64, ktop: Int64, kbot: Int64, nshfts: Int64,
        sr: Array<Float64>, si: Array<Float64>, h: Array<Float64>, ldh: Int64, iloz: Int64, ihiz: Int64,
        z: Array<Float64>, ldz: Int64, v: Array<Float64>, ldv: Int64, u: Array<Float64>, ldu: Int64, nv: Int64,
        wv: Array<Float64>, ldwv: Int64, nh: Int64, wh: Array<Float64>, ldwh: Int64): Unit {
        if (kacc22 != 0 && kacc22 != 1 && kacc22 != 2) {
            throw Exception(badKacc22)
        } else if (n < 0) {
            throw Exception(nLT0)
        } else if (ktop < 0 || n <= ktop) {
            throw Exception(badKtop)
        } else if (kbot < 0 || n <= kbot) {
            throw Exception(badKbot)
        } else if (nshfts < 0) {
            throw Exception(nshftsLT0)
        } else if (nshfts % 2 != 0) {
            throw Exception(nshftsOdd)
        } else if (sr.size != nshfts) {
            throw Exception(badLenSr)
        } else if (si.size != nshfts) {
            throw Exception(badLenSi)
        } else if (ldh < max(1, n)) {
            throw Exception(badLdH)
        } else if (h.size < (n - 1) * ldh + n) {
            throw Exception(shortH)
        } else if (wantz && ihiz >= n) {
            throw Exception(badIhiz)
        } else if (wantz && iloz < 0 || ihiz < iloz) {
            throw Exception(badIloz)
        } else if (ldz < 1 || wantz && ldz < n) {
            throw Exception(badLdZ)
        } else if (wantz && z.size < (n - 1) * ldz + n) {
            throw Exception(shortZ)
        } else if (ldv < 3) {
            throw Exception(badLdV)
        } else if (v.size < (nshfts / 2 - 1) * ldv + 3) {
            throw Exception(shortV)
        } else if (ldu < max(1, 2 * nshfts)) {
            throw Exception(badLdU)
        } else if (u.size < (2 * nshfts - 1) * ldu + 2 * nshfts) {
            throw Exception(shortU)
        } else if (nv < 0) {
            throw Exception(nvLT0)
        } else if (ldwv < max(1, 2 * nshfts)) {
            throw Exception(badLdWV)
        } else if (wv.size < (nv - 1) * ldwv + 2 * nshfts) {
            throw Exception(shortWV)
        } else if (nh < 0) {
            throw Exception(nhLT0)
        } else if (ldwh < max(1, nh)) {
            throw Exception(badLdWH)
        } else if (wh.size < (2 * nshfts - 1) * ldwh + nh) {
            throw Exception(shortWH)
        } else if (ktop > 0 && h[ktop * ldh + ktop - 1] != 0.0) {
            throw Exception(notIsolated)
        } else if (kbot < n - 1 && h[(kbot + 1) * ldh + kbot] != 0.0) {
            throw Exception(notIsolated)
        }

        if (nshfts < 2) {
            return
        }

        if (ktop >= kbot) {
            return
        }

        for (i in 0..nshfts - 2 : 2) {
            if (si[i] == -si[i + 1]) {
                continue
            }
            (sr[i], sr[i + 1], sr[i + 2]) = (sr[i + 1], sr[i + 2], sr[i])
            (si[i], si[i + 1], si[i + 2]) = (si[i + 1], si[i + 2], si[i])
        }

        let ns = nshfts

        let safmin = dlamchS
        let ulp = dlamchP
        let smlnum = safmin * Float64(n) / ulp

        let accum = kacc22 == 1 || kacc22 == 2

        if (ktop + 2 <= kbot) {
            h[(ktop + 2) * ldh + ktop] = 0.0
        }

        let nbmps = ns / 2

        let kdu = 4 * nbmps

        for (incol in ktop - 2 * nbmps + 1..=kbot - 2 : 2 * nbmps) {
            var jtop: Int64 = 0
            if (accum) {
                jtop = max(ktop, incol)
            } else if (wantt) {
            } else {
                jtop = ktop
            }
            let ndcol = incol + kdu
            if (accum) {
                dlaset(All, kdu, kdu, 0.0, 1.0, u, ldu)
            }

            for (krcol in incol..=min(incol + 2 * nbmps - 1, kbot - 2)) {
                let mtop = max(0, (ktop - krcol) / 2)
                let mbot = min(nbmps, (kbot - krcol - 1) / 2) - 1
                let m22 = mbot + 1
                let bmp22 = (mbot < nbmps - 1) && (krcol + 2 * m22 == kbot - 2)

                if (bmp22) {
                    let k = krcol + 2 * m22
                    if (k == ktop - 1) {
                        dlaqr1(2, h[(k + 1) * ldh + k + 1..], ldh, sr[2 * m22], si[2 * m22], sr[2 * m22 + 1],
                            si[2 * m22 + 1], v[m22 * ldv..m22 * ldv + 2])
                        let beta = v[m22 * ldv]
                        (_, v[m22 * ldv]) = dlarfg(2, beta, v[m22 * ldv + 1..m22 * ldv + 2], 1)
                    } else {
                        var beta = h[(k + 1) * ldh + k]
                        v[m22 * ldv + 1] = h[(k + 2) * ldh + k]
                        (beta, v[m22 * ldv]) = dlarfg(2, beta, v[m22 * ldv + 1..m22 * ldv + 2], 1)
                        h[(k + 1) * ldh + k] = beta
                        h[(k + 2) * ldh + k] = 0.0
                    }

                    var t1 = v[m22 * ldv]
                    var t2 = t1 * v[m22 * ldv + 1]
                    for (j in jtop..=min(kbot, k + 3)) {
                        let refsum = h[j * ldh + k + 1] + v[m22 * ldv + 1] * h[j * ldh + k + 2]
                        h[j * ldh + k + 1] -= refsum * t1
                        h[j * ldh + k + 2] -= refsum * t2
                    }

                    var jbot: Int64
                    if (accum) {
                        jbot = min(ndcol, kbot)
                    } else if (wantt) {
                        jbot = n - 1
                    } else {
                        jbot = kbot
                    }
                    t1 = v[m22 * ldv]
                    t2 = t1 * v[m22 * ldv + 1]
                    for (j in k + 1..=jbot) {
                        let refsum = h[(k + 1) * ldh + j] + v[m22 * ldv + 1] * h[(k + 2) * ldh + j]
                        h[(k + 1) * ldh + j] -= refsum * t1
                        h[(k + 2) * ldh + j] -= refsum * t2
                    }

                    if (k >= ktop && h[(k + 1) * ldh + k] != 0.0) {
                        var tst1 = abs(h[k * ldh + k]) + abs(h[(k + 1) * ldh + k + 1])
                        if (tst1 == 0.0) {
                            if (k >= ktop + 1) {
                                tst1 += abs(h[k * ldh + k - 1])
                            }
                            if (k >= ktop + 2) {
                                tst1 += abs(h[k * ldh + k - 2])
                            }
                            if (k >= ktop + 3) {
                                tst1 += abs(h[k * ldh + k - 3])
                            }
                            if (k <= kbot - 2) {
                                tst1 += abs(h[(k + 2) * ldh + k + 1])
                            }
                            if (k <= kbot - 3) {
                                tst1 += abs(h[(k + 3) * ldh + k + 1])
                            }
                            if (k <= kbot - 4) {
                                tst1 += abs(h[(k + 4) * ldh + k + 1])
                            }
                        }
                        if (abs(h[(k + 1) * ldh + k]) <= max(smlnum, ulp * tst1)) {
                            let h12 = max(abs(h[(k + 1) * ldh + k]), abs(h[k * ldh + k + 1]))
                            let h21 = min(abs(h[(k + 1) * ldh + k]), abs(h[k * ldh + k + 1]))
                            let h11 = max(abs(h[(k + 1) * ldh + k + 1]), abs(h[k * ldh + k] - h[(k + 1) * ldh + k + 1]))
                            let h22 = min(abs(h[(k + 1) * ldh + k + 1]), abs(h[k * ldh + k] - h[(k + 1) * ldh + k + 1]))
                            let scl = h11 + h12
                            let tst2 = h22 * (h11 / scl)
                            if (tst2 == 0.0 || h21 * (h12 / scl) <= max(smlnum, ulp * tst2)) {
                                h[(k + 1) * ldh + k] = 0.0
                            }
                        }
                    }

                    if (accum) {
                        let kms = k - incol - 1
                        let t1 = v[m22 * ldv]
                        let t2 = t1 * v[m22 * ldv + 1]
                        for (j in max(0, ktop - incol - 1)..kdu - 1) {
                            let refsum = u[j * ldu + kms + 1] + v[m22 * ldv + 1] * u[j * ldu + kms + 2]
                            u[j * ldu + kms + 1] -= refsum * t1
                            u[j * ldu + kms + 2] -= refsum * t2
                        }
                    } else if (wantz) {
                        let t1 = v[m22 * ldv]
                        let t2 = t1 * v[m22 * ldv + 1]
                        for (j in iloz..=ihiz) {
                            let refsum = z[j * ldz + k + 1] + v[m22 * ldv + 1] * z[j * ldz + k + 2]
                            z[j * ldz + k + 1] -= refsum * t1
                            z[j * ldz + k + 2] -= refsum * t2
                        }
                    }
                }

                for (m in mbot..=mtop : -1) {
                    let k = krcol + 2 * m
                    if (k == ktop - 1) {
                        dlaqr1(3, h[ktop * ldh + ktop..], ldh, sr[2 * m], si[2 * m], sr[2 * m + 1], si[2 * m + 1],
                            v[m * ldv..m * ldv + 3])
                        let alpha = v[m * ldv]
                        (_, v[m * ldv]) = dlarfg(3, alpha, v[m * ldv + 1..m * ldv + 3], 1)
                    } else {
                        var t1 = v[m * ldv]
                        var t2 = t1 * v[m * ldv + 1]
                        var t3 = t1 * v[m * ldv + 2]
                        var refsum = v[m * ldv + 2] * h[(k + 3) * ldh + k + 2]
                        h[(k + 3) * ldh + k] = -refsum * t1
                        h[(k + 3) * ldh + k + 1] = -refsum * t2
                        h[(k + 3) * ldh + k + 2] -= refsum * t3

                        var beta = h[(k + 1) * ldh + k]
                        v[m * ldv + 1] = h[(k + 2) * ldh + k]
                        v[m * ldv + 2] = h[(k + 3) * ldh + k]
                        (beta, v[m * ldv]) = dlarfg(3, beta, v[m * ldv + 1..m * ldv + 3], 1)

                        if (h[(k + 3) * ldh + k] != 0.0 || h[(k + 3) * ldh + k + 1] != 0.0 || h[(k + 3) * ldh + k + 2] ==
                            0.0) {
                            h[(k + 1) * ldh + k] = beta
                            h[(k + 2) * ldh + k] = 0.0
                            h[(k + 3) * ldh + k] = 0.0
                        } else {
                            var vt = Array<Float64>(3, item: 0.0)
                            dlaqr1(3, h[(k + 1) * ldh + k + 1..], ldh, sr[2 * m], si[2 * m], sr[2 * m + 1], si[2 * m + 1],
                                vt[..])
                            (_, vt[0]) = dlarfg(3, vt[0], vt[1..3], 1)
                            t1 = vt[0]
                            t2 = t1 * vt[1]
                            t3 = t1 * vt[2]
                            refsum = h[(k + 1) * ldh + k] + vt[1] * h[(k + 2) * ldh + k]
                            let dsum = abs(h[k * ldh + k]) + abs(h[(k + 1) * ldh + k + 1]) + abs(h[(k + 2) * ldh + k + 2]
                            )
                            if (abs(h[(k + 2) * ldh + k] - refsum * t2) + abs(refsum * t3) > ulp * dsum) {
                                h[(k + 1) * ldh + k] = beta
                                h[(k + 2) * ldh + k] = 0.0
                                h[(k + 3) * ldh + k] = 0.0
                            } else {
                                h[(k + 1) * ldh + k] -= refsum * t1
                                h[(k + 2) * ldh + k] = 0.0
                                h[(k + 3) * ldh + k] = 0.0
                                v[m * ldv] = vt[0]
                                v[m * ldv + 1] = vt[1]
                                v[m * ldv + 2] = vt[2]
                            }
                        }
                    }

                    let t1 = v[m * ldv]
                    let t2 = t1 * v[m * ldv + 1]
                    let t3 = t1 * v[m * ldv + 2]
                    for (j in jtop..=min(kbot, k + 3)) {
                        let refsum = h[j * ldh + k + 1] + v[m * ldv + 1] * h[j * ldh + k + 2] + v[m * ldv + 2] * h[j *
                                ldh + k + 3]
                        h[j * ldh + k + 1] -= refsum * t1
                        h[j * ldh + k + 2] -= refsum * t2
                        h[j * ldh + k + 3] -= refsum * t3
                    }

                    let refsum = h[(k + 1) * ldh + k + 1] + v[m * ldv + 1] * h[(k + 2) * ldh + k + 1] + v[m * ldv + 2] *
                        h[(k + 3) * ldh + k + 1]
                    h[(k + 1) * ldh + k + 1] -= refsum * t1
                    h[(k + 2) * ldh + k + 1] -= refsum * t2
                    h[(k + 3) * ldh + k + 1] -= refsum * t3

                    if (k < ktop) {
                        continue
                    }
                    if (h[(k + 1) * ldh + k] != 0.0) {
                        var tst1 = abs(h[k * ldh + k]) + abs(h[(k + 1) * ldh + k + 1])
                        if (tst1 == 0.0) {
                            if (k >= ktop + 1) {
                                tst1 += abs(h[k * ldh + k - 1])
                            }
                            if (k >= ktop + 2) {
                                tst1 += abs(h[k * ldh + k - 2])
                            }
                            if (k >= ktop + 3) {
                                tst1 += abs(h[k * ldh + k - 3])
                            }
                            if (k <= kbot - 2) {
                                tst1 += abs(h[(k + 2) * ldh + k + 1])
                            }
                            if (k <= kbot - 3) {
                                tst1 += abs(h[(k + 3) * ldh + k + 1])
                            }
                            if (k <= kbot - 4) {
                                tst1 += abs(h[(k + 4) * ldh + k + 1])
                            }
                        }
                        if (abs(h[(k + 1) * ldh + k]) <= max(smlnum, ulp * tst1)) {
                            let h12 = max(abs(h[(k + 1) * ldh + k]), abs(h[k * ldh + k + 1]))
                            let h21 = min(abs(h[(k + 1) * ldh + k]), abs(h[k * ldh + k + 1]))
                            let h11 = max(abs(h[(k + 1) * ldh + k + 1]), abs(h[k * ldh + k] - h[(k + 1) * ldh + k + 1]))
                            let h22 = min(abs(h[(k + 1) * ldh + k + 1]), abs(h[k * ldh + k] - h[(k + 1) * ldh + k + 1]))
                            let scl = h11 + h12
                            let tst2 = h22 * (h11 / scl)
                            if (tst2 == 0.0 || h21 * (h12 / scl) <= max(smlnum, ulp * tst2)) {
                                h[(k + 1) * ldh + k] = 0.0
                            }
                        }
                    }
                }

                var jbot: Int64
                if (accum) {
                    jbot = min(ndcol, kbot)
                } else if (wantt) {
                    jbot = n - 1
                } else {
                    jbot = kbot
                }
                for (m in mbot..=mtop : -1) {
                    let k = krcol + 2 * m
                    let t1 = v[m * ldv]
                    let t2 = t1 * v[m * ldv + 1]
                    let t3 = t1 * v[m * ldv + 2]
                    for (j in max(ktop, krcol + 2 * (m + 1))..=jbot) {
                        let refsum = h[(k + 1) * ldh + j] + v[m * ldv + 1] * h[(k + 2) * ldh + j] + v[m * ldv + 2] * h[(k +
                                3) * ldh + j]
                        h[(k + 1) * ldh + j] -= refsum * t1
                        h[(k + 2) * ldh + j] -= refsum * t2
                        h[(k + 3) * ldh + j] -= refsum * t3
                    }
                }

                if (accum) {
                    for (m in mbot..=mtop : -1) {
                        let k = krcol + 2 * m
                        let kms = k - incol - 1
                        var i2 = max(0, ktop - incol - 1)
                        i2 = max(i2, kms - (krcol - incol))
                        let i4 = min(kdu, krcol + 2 * mbot - incol + 5)
                        let t1 = v[m * ldv]
                        let t2 = t1 * v[m * ldv + 1]
                        let t3 = t1 * v[m * ldv + 2]
                        for (j in i2..i4 - 1) {
                            let refsum = u[j * ldu + kms + 1] + v[m * ldv + 1] * u[j * ldu + kms + 2] + v[m * ldv + 2] *
                                u[j * ldu + kms + 3]
                            u[j * ldu + kms + 1] -= refsum * t1
                            u[j * ldu + kms + 2] -= refsum * t2
                            u[j * ldu + kms + 3] -= refsum * t3
                        }
                    }
                } else if (wantz) {
                    for (m in mbot..=mtop : -1) {
                        let k = krcol + 2 * m
                        let t1 = v[m * ldv]
                        let t2 = t1 * v[m * ldv + 1]
                        let t3 = t1 * v[m * ldv + 2]
                        for (j in iloz..=ihiz) {
                            let refsum = z[j * ldz + k + 1] + v[m * ldv + 1] * z[j * ldz + k + 2] + v[m * ldv + 2] * z[j *
                                    ldz + k + 3]
                            z[j * ldz + k + 1] -= refsum * t1
                            z[j * ldz + k + 2] -= refsum * t2
                            z[j * ldz + k + 3] -= refsum * t3
                        }
                    }
                }
            }

            if (!accum) {
                continue
            }
            jtop = ktop
            var jbot = kbot
            if (wantt) {
                jtop = 0
                jbot = n - 1
            }
            let nFloat64 = nFloat64Implementation()
            let k1 = max(0, ktop - incol - 1)
            let nu = kdu - max(0, ndcol - kbot) - k1

            for (jcol in min(ndcol, kbot) + 1..=jbot : nh) {
                let jlen = min(nh, jbot - jcol + 1)
                nFloat64.dgemm(Trans, NoTrans, nu, jlen, nu, 1.0, u[k1 * ldu + k1..], ldu,
                    h[(incol + k1 + 1) * ldh + jcol..], ldh, 0.0, wh, ldwh)
                dlacpy(All, nu, jlen, wh, ldwh, h[(incol + k1 + 1) * ldh + jcol..], ldh)
            }

            for (jrow in jtop..max(ktop, incol) : nv) {
                let jlen = min(nv, max(ktop, incol) - jrow)
                nFloat64.dgemm(NoTrans, NoTrans, jlen, nu, nu, 1.0, h[jrow * ldh + incol + k1 + 1..], ldh,
                    u[k1 * ldu + k1..], ldu, 0.0, wv, ldwv)
                dlacpy(All, jlen, nu, wv, ldwv, h[jrow * ldh + incol + k1 + 1..], ldh)
            }

            if (wantz) {
                for (jrow in iloz..=ihiz : nv) {
                    let jlen = min(nv, ihiz - jrow + 1)
                    nFloat64.dgemm(NoTrans, NoTrans, jlen, nu, nu, 1.0, z[jrow * ldz + incol + k1 + 1..], ldz,
                        u[k1 * ldu + k1..], ldu, 0.0, wv, ldwv)
                    dlacpy(All, jlen, nu, wv, ldwv, z[jrow * ldz + incol + k1 + 1..], ldz)
                }
            }
        }
    }
}
