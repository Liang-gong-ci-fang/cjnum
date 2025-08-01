package cjnum.lapack.cjnum

import cjnum.blas.blas64.*
import cjnum.blas.*
import cjnum.lapack.*
import cjnum.util.*
import std.math.*

extend LFloat64Implementation {
    public func dtrevc3(side: EVSide, howmny: EVHowMany, selected: Array<Bool>, n: Int64, t: Array<Float64>, ldt: Int64,
        vl: Array<Float64>, ldvl: Int64, vr: Array<Float64>, ldvr: Int64, mm: Int64, work: Array<Float64>, lwork: Int64): Int64 {
        var m: Int64 = 0
        // let bothv = side == EVBoth
        let rightv = side == EVRight || (side == EVBoth)
        let leftv = side == EVLeft || (side == EVBoth)
        if (!rightv && !leftv) {
            throw Exception(badEVSide)
        } else if (howmny != EVAll && howmny != EVAllMulQ && howmny != EVSelected) {
            throw Exception(badEVHowMany)
        } else if (n < 0) {
            throw Exception(nLT0)
        } else if (ldt < max(1, n)) {
            throw Exception(badLdT)
        } else if (mm < 0) {
            throw Exception(mmLT0)
        } else if (ldvl < 1) {
            throw Exception(badLdVL)
        } else if (ldvr < 1) {
            throw Exception(badLdVR)
        } else if (lwork < max(1, 3 * n) && lwork != -1) {
            throw Exception(badLWork)
        } else if (work.size < max(1, lwork)) {
            throw Exception(shortWork)
        }

        if (n == 0) {
            work[0] = 1.0
            return 0
        }

        if (t.size < (n - 1) * ldt + n) {
            throw Exception(shortT)
        }

        if (howmny == EVSelected) {
            if (selected.size != n) {
                throw Exception(badLenSelected)
            }
            var j = 0
            while (j < n) {
                if (j == n - 1 || t[(j + 1) * ldt + j] == 0.0) {
                    if (selected[j]) {
                        m++
                    }
                    j++
                } else {
                    if (selected[j] || selected[j + 1]) {
                        selected[j] = true
                        selected[j + 1] = false
                        m += 2
                    }
                    j += 2
                }
            }
        } else {
            m = n
        }
        if (mm < m) {
            throw Exception(badMm)
        }

        var nb = ilaenv(1, "DTREVC", side.toString() + howmny.toString(), n, -1, -1, -1)
        if (lwork == -1) {
            work[0] = Float64(n + 2 * n * nb)
            return m
        }

        if (m == 0) {
            return 0
        }

        if (leftv && ldvl < mm) {
            throw Exception(badLdVL)
        } else if (leftv && vl.size < (n - 1) * ldvl + mm) {
            throw Exception(shortVL)
        } else if (rightv && ldvr < mm) {
            throw Exception(badLdVR)
        } else if (rightv && vr.size < (n - 1) * ldvr + mm) {
            throw Exception(shortVR)
        }

        const nbmin = 8
        const nbmax = 128

        if (howmny == EVAllMulQ && lwork >= n + 2 * n * nbmin) {
            nb = min((lwork - n) / (2 * n), nbmax)
            dlaset(All, n, 1 + 2 * nb, 0.0, 0.0, work[..n + 2 * nb * n], 1 + 2 * nb)
        } else {
            nb = 1
        }

        let ulp = dlamchP
        let smlnum = Float64(n) / ulp * dlamchS
        let bignum = (1.0 - ulp) / smlnum

        let norms = work[..n]
        let ldb = 2 * nb
        let b = work[n..n + n * ldb]

        norms[0] = 0.0
        for (j in 1..n) {
            var cn: Float64 = 0.0
            for (i in 0..j) {
                cn += abs(t[i * ldt + j])
            }
            norms[j] = cn
        }

        let nFloat64 = nFloat64Implementation()

        let x: Array<Float64> = Array<Float64>(4, repeat: 0.0)
        var iv: Int64
        var is_: Int64
        var ip: Int64
        let iscomplex: Array<Int64> = Array<Int64>(nbmax, repeat: 0)

        var leftev = true
        if (side == EVLeft) {
            leftev = false
        }

        if (leftev) {
            iv = max(2, nb) - 1
            ip = 0
            is_ = m - 1
            for (ki in n - 1..=0 : -1) {
                if (ip == -1) {
                    ip = 1
                    continue
                }

                if (ki == 0 || t[ki * ldt + ki - 1] == 0.0) {
                    ip = 0
                } else {
                    ip = -1
                }

                if (howmny == EVSelected) {
                    if (ip == 0) {
                        if (!selected[ki]) {
                            continue
                        }
                    } else if (!selected[ki - 1]) {
                        continue
                    }
                }

                let wr = t[ki * ldt + ki]
                var wi: Float64 = 0.0
                if (ip != 0) {
                    wi = sqrt(abs(t[ki * ldt + ki - 1])) * sqrt(abs(t[(ki - 1) * ldt + ki]))
                }
                let smin = nanMax(ulp * (abs(wr) + abs(wi)), smlnum)

                if (ip == 0) {
                    b[ki * ldb + iv] = 1.0
                    for (k in 0..ki) {
                        b[k * ldb + iv] = -t[k * ldt + ki]
                    }
                    var j = ki - 1
                    while (j >= 0) {
                        if (j == 0 || t[j * ldt + j - 1] == 0.0) {
                            var (scale, xnorm, _) = dlaln2(false, 1, 1, smin, 1.0, t[j * ldt + j..], ldt, 1.0, 1.0,
                                b[j * ldb + iv..], ldb, wr, 0.0, x[..1], 2)
                            if (xnorm > 1.0 && norms[j] > bignum / xnorm) {
                                x[0] /= xnorm
                                scale /= xnorm
                            }
                            if (scale != 1.0) {
                                nFloat64.dscal(ki + 1, scale, b[iv..], ldb)
                            }
                            b[j * ldb + iv] = x[0]
                            nFloat64.daxpy(j, -x[0], t[j..], ldt, b[iv..], ldb)
                            j--
                        } else {
                            var (scale, xnorm, _) = dlaln2(false, 2, 1, smin, 1.0, t[(j - 1) * ldt + j - 1..], ldt, 1.0,
                                1.0, b[(j - 1) * ldb + iv..], ldb, wr, 0.0, x[..3], 2)
                            if (xnorm > 1.0) {
                                let beta = nanMax(norms[j - 1], norms[j])
                                if (beta > bignum / xnorm) {
                                    x[0] /= xnorm
                                    x[2] /= xnorm
                                    scale /= xnorm
                                }
                            }
                            if (scale != 1.0) {
                                nFloat64.dscal(ki + 1, scale, b[iv..], ldb)
                            }
                            b[(j - 1) * ldb + iv] = x[0]
                            b[j * ldb + iv] = x[2]
                            nFloat64.daxpy(j - 1, -x[0], t[j - 1..], ldt, b[iv..], ldb)
                            nFloat64.daxpy(j - 1, -x[2], t[j..], ldt, b[iv..], ldb)
                            j -= 2
                        }
                    }
                    if (howmny != EVAllMulQ) {
                        nFloat64.dcopy(ki + 1, b[iv..], ldb, vr[is_..], ldvr)
                        let ii = nFloat64.idamax(ki + 1, vr[is_..], ldvr)
                        let remax = 1.0 / abs(vr[ii * ldvr + is_])
                        nFloat64.dscal(ki + 1, remax, vr[is_..], ldvr)
                        for (k in ki + 1..n) {
                            vr[k * ldvr + is_] = 0.0
                        }
                    } else if (nb == 1) {
                        if (ki > 0) {
                            nFloat64.dgemv(NoTrans, n, ki, 1.0, vr, ldvr, b[iv..], ldb, b[ki * ldb + iv], vr[ki..], ldvr
                            )
                        }
                        let ii = nFloat64.idamax(n, vr[ki..], ldvr)
                        let remax = 1.0 / abs(vr[ii * ldvr + ki])
                        nFloat64.dscal(n, remax, vr[ki..], ldvr)
                    } else {
                        for (k in ki + 1..n) {
                            b[k * ldb + iv] = 0.0
                        }
                        iscomplex[iv] = ip
                    }
                } else {
                    if (abs(t[(ki - 1) * ldt + ki]) >= abs(t[ki * ldt + ki - 1])) {
                        b[(ki - 1) * ldb + iv - 1] = 1.0
                        b[ki * ldb + iv] = wi / t[(ki - 1) * ldt + ki]
                    } else {
                        b[(ki - 1) * ldb + iv - 1] = -wi / t[ki * ldt + ki - 1]
                        b[ki * ldb + iv] = 1.0
                    }
                    b[ki * ldb + iv - 1] = 0.0
                    b[(ki - 1) * ldb + iv] = 0.0
                    for (k in 0..ki - 1) {
                        b[k * ldb + iv - 1] = -b[(ki - 1) * ldb + iv - 1] * t[k * ldt + ki - 1]
                        b[k * ldb + iv] = -b[ki * ldb + iv] * t[k * ldt + ki]
                    }
                    var j = ki - 2
                    while (j >= 0) {
                        if (j == 0 || t[j * ldt + j - 1] == 0.0) {
                            var (scale, xnorm, _) = dlaln2(false, 1, 2, smin, 1.0, t[j * ldt + j..], ldt, 1.0, 1.0,
                                b[j * ldb + iv - 1..], ldb, wr, wi, x[..2], 2)
                            if (xnorm > 1.0 && norms[j] > bignum / xnorm) {
                                x[0] /= xnorm
                                x[1] /= xnorm
                                scale /= xnorm
                            }
                            if (scale != 1.0) {
                                nFloat64.dscal(ki + 1, scale, b[iv - 1..], ldb)
                                nFloat64.dscal(ki + 1, scale, b[iv..], ldb)
                            }
                            b[j * ldb + iv - 1] = x[0]
                            b[j * ldb + iv] = x[1]
                            nFloat64.daxpy(j, -x[0], t[j..], ldt, b[iv - 1..], ldb)
                            nFloat64.daxpy(j, -x[1], t[j..], ldt, b[iv..], ldb)
                            j--
                        } else {
                            var (scale, xnorm, _) = dlaln2(false, 2, 2, smin, 1.0, t[(j - 1) * ldt + j - 1..], ldt, 1.0,
                                1.0, b[(j - 1) * ldb + iv - 1..], ldb, wr, wi, x[..], 2)
                            if (xnorm > 1.0) {
                                let beta = nanMax(norms[j - 1], norms[j])
                                if (beta > bignum / xnorm) {
                                    let rec = 1.0 / xnorm
                                    x[0] *= rec
                                    x[1] *= rec
                                    x[2] *= rec
                                    x[3] *= rec
                                    scale *= rec
                                }
                            }
                            if (scale != 1.0) {
                                nFloat64.dscal(ki + 1, scale, b[iv - 1..], ldb)
                                nFloat64.dscal(ki + 1, scale, b[iv..], ldb)
                            }
                            b[(j - 1) * ldb + iv - 1] = x[0]
                            b[(j - 1) * ldb + iv] = x[1]
                            b[j * ldb + iv - 1] = x[2]
                            b[j * ldb + iv] = x[3]
                            nFloat64.daxpy(j - 1, -x[0], t[j - 1..], ldt, b[iv - 1..], ldb)
                            nFloat64.daxpy(j - 1, -x[1], t[j - 1..], ldt, b[iv..], ldb)
                            nFloat64.daxpy(j - 1, -x[2], t[j..], ldt, b[iv - 1..], ldb)
                            nFloat64.daxpy(j - 1, -x[3], t[j..], ldt, b[iv..], ldb)
                            j -= 2
                        }
                    }
                    if (howmny != EVAllMulQ) {
                        nFloat64.dcopy(ki + 1, b[iv - 1..], ldb, vr[is_ - 1..], ldvr)
                        nFloat64.dcopy(ki + 1, b[iv..], ldb, vr[is_..], ldvr)
                        var emax = 0.0
                        for (k in 0..=ki) {
                            emax = nanMax(emax, abs(vr[k * ldvr + is_ - 1]) + abs(vr[k * ldvr + is_]))
                        }
                        let remax = 1.0 / emax
                        nFloat64.dscal(ki + 1, remax, vr[is_ - 1..], ldvr)
                        nFloat64.dscal(ki + 1, remax, vr[is_..], ldvr)
                        for (k in ki + 1..n) {
                            vr[k * ldvr + is_ - 1] = 0.0
                            vr[k * ldvr + is_] = 0.0
                        }
                    } else if (nb == 1) {
                        if (ki - 1 > 0) {
                            nFloat64.dgemv(NoTrans, n, ki - 1, 1.0, vr, ldvr, b[iv - 1..], ldb, b[(ki - 1) * ldb + iv - 1],
                                vr[ki - 1..], ldvr)
                            nFloat64.dgemv(NoTrans, n, ki - 1, 1.0, vr, ldvr, b[iv..], ldb, b[ki * ldb + iv], vr[ki..],
                                ldvr)
                        } else {
                            nFloat64.dscal(n, b[(ki - 1) * ldb + iv - 1], vr[ki - 1..], ldvr)
                            nFloat64.dscal(n, b[ki * ldb + iv], vr[ki..], ldvr)
                        }
                        var emax = 0.0
                        for (k in 0..n) {
                            emax = nanMax(emax, abs(vr[k * ldvr + ki - 1]) + abs(vr[k * ldvr + ki]))
                        }
                        let remax = 1.0 / emax
                        nFloat64.dscal(n, remax, vr[ki - 1..], ldvr)
                        nFloat64.dscal(n, remax, vr[ki..], ldvr)
                    } else {
                        for (k in ki + 1..n) {
                            b[k * ldb + iv - 1] = 0.0
                            b[k * ldb + iv] = 0.0
                        }
                        iscomplex[iv - 1] = -ip
                        iscomplex[iv] = ip
                        iv--
                    }
                }
                if (nb > 1) {
                    var ki2 = ki
                    if (ip != 0) {
                        ki2--
                    }
                    if (iv < 2 || ki2 == 0) {
                        nFloat64.dgemm(NoTrans, NoTrans, n, nb - iv, ki2 + nb - iv, 1.0, vr, ldvr, b[iv..], ldb, 0.0,
                            b[nb + iv..], ldb)
                        var remax: Float64 = 0.0
                        for (k in iv..nb) {
                            if (iscomplex[k] == 0) {
                                let ii = nFloat64.idamax(n, b[nb + k..], ldb)
                                remax = 1.0 / abs(b[ii * ldb + nb + k])
                            } else if (iscomplex[k] == 1) {
                                var emax = 0.0
                                for (ii in 0..n) {
                                    emax = nanMax(emax, abs(b[ii * ldb + nb + k]) + abs(b[ii * ldb + nb + k + 1]))
                                }
                                remax = 1.0 / emax
                            }
                            nFloat64.dscal(n, remax, b[nb + k..], ldb)
                        }
                        dlacpy(All, n, nb - iv, b[nb + iv..], ldb, vr[ki2..], ldvr)
                        iv = nb - 1
                    } else {
                        iv--
                    }
                }
                is_--
                if (ip != 0) {
                    is_--
                }
            }

            if (side == EVRight) {
                return m
            }
        }

        iv = 0
        ip = 0
        is_ = 0
        for (ki in 0..n) {
            if (ip == 1) {
                ip = -1
                continue
            }

            if (ki == n - 1 || t[(ki + 1) * ldt + ki] == 0.0) {
                ip = 0
            } else {
                ip = 1
            }
            if (howmny == EVSelected && !selected[ki]) {
                continue
            }

            let wr = t[ki * ldt + ki]
            var wi: Float64 = 0.0
            if (ip != 0) {
                wi = sqrt(abs(t[ki * ldt + ki + 1])) * sqrt(abs(t[(ki + 1) * ldt + ki]))
            }
            let smin = nanMax(ulp * (abs(wr) + abs(wi)), smlnum)

            if (ip == 0) {
                b[ki * ldb + iv] = 1.0
                for (k in ki + 1..n) {
                    b[k * ldb + iv] = -t[ki * ldt + k]
                }
                var vmax = 1.0
                var vcrit = bignum
                var j = ki + 1
                while (j < n) {
                    if (j == n - 1 || t[(j + 1) * ldt + j] == 0.0) {
                        if (norms[j] > vcrit) {
                            let rec = 1.0 / vmax
                            nFloat64.dscal(n - ki, rec, b[ki * ldb + iv..], ldb)
                            vmax = 1.0
                        }
                        b[j * ldb + iv] -= nFloat64.ddot(j - ki - 1, t[(ki + 1) * ldt + j..], ldt, b[(ki + 1) * ldb + iv..],
                            ldb)
                        let (scale, _, _) = dlaln2(false, 1, 1, smin, 1.0, t[j * ldt + j..], ldt, 1.0, 1.0,
                            b[j * ldb + iv..], ldb, wr, 0.0, x[..1], 2)
                        if (scale != 1.0) {
                            nFloat64.dscal(n - ki, scale, b[ki * ldb + iv..], ldb)
                        }
                        b[j * ldb + iv] = x[0]
                        vmax = nanMax(abs(b[j * ldb + iv]), vmax)
                        vcrit = bignum / vmax
                        j++
                    } else {
                        let beta = nanMax(norms[j], norms[j + 1])
                        if (beta > vcrit) {
                            nFloat64.dscal(n - ki, 1.0 / vmax, b[ki * ldb + iv..], ldb)
                            vmax = 1.0
                        }
                        b[j * ldb + iv] -= nFloat64.ddot(j - ki - 1, t[(ki + 1) * ldt + j..], ldt, b[(ki + 1) * ldb + iv..],
                            ldb)
                        b[(j + 1) * ldb + iv] -= nFloat64.ddot(j - ki - 1, t[(ki + 1) * ldt + j + 1..], ldt,
                            b[(ki + 1) * ldb + iv..], ldb)
                        let (scale, _, _) = dlaln2(true, 2, 1, smin, 1.0, t[j * ldt + j..], ldt, 1.0, 1.0,
                            b[j * ldb + iv..], ldb, wr, 0.0, x[..3], 2)
                        if (scale != 1.0) {
                            nFloat64.dscal(n - ki, scale, b[ki * ldb + iv..], ldb)
                        }
                        b[j * ldb + iv] = x[0]
                        b[(j + 1) * ldb + iv] = x[2]
                        vmax = nanMax(vmax, nanMax(abs(b[j * ldb + iv]), abs(b[(j + 1) * ldb + iv])))
                        vcrit = bignum / vmax
                        j += 2
                    }
                }
                if (howmny != EVAllMulQ) {
                    nFloat64.dcopy(n - ki, b[ki * ldb + iv..], ldb, vl[ki * ldvl + is_..], ldvl)
                    let ii = nFloat64.idamax(n - ki, vl[ki * ldvl + is_..], ldvl) + ki
                    let remax = 1.0 / abs(vl[ii * ldvl + is_])
                    nFloat64.dscal(n - ki, remax, vl[ki * ldvl + is_..], ldvl)
                    for (k in 0..ki) {
                        vl[k * ldvl + is_] = 0.0
                    }
                } else if (nb == 1) {
                    if (n - ki - 1 > 0) {
                        nFloat64.dgemv(NoTrans, n, n - ki - 1, 1.0, vl[ki + 1..], ldvl, b[(ki + 1) * ldb + iv..], ldb,
                            b[ki * ldb + iv], vl[ki..], ldvl)
                    }
                    let ii = nFloat64.idamax(n, vl[ki..], ldvl)
                    let remax = 1.0 / abs(vl[ii * ldvl + ki])
                    nFloat64.dscal(n, remax, vl[ki..], ldvl)
                } else {
                    for (k in 0..ki) {
                        b[k * ldb + iv] = 0.0
                    }
                    iscomplex[iv] = ip
                }
            } else {
                if (abs(t[ki * ldt + ki + 1]) >= abs(t[(ki + 1) * ldt + ki])) {
                    b[ki * ldb + iv] = wi / t[ki * ldt + ki + 1]
                    b[(ki + 1) * ldb + iv + 1] = 1.0
                } else {
                    b[ki * ldb + iv] = 1.0
                    b[(ki + 1) * ldb + iv + 1] = -wi / t[(ki + 1) * ldt + ki]
                }
                b[(ki + 1) * ldb + iv] = 0.0
                b[ki * ldb + iv + 1] = 0.0
                for (k in ki + 2..n) {
                    b[k * ldb + iv] = -b[ki * ldb + iv] * t[ki * ldt + k]
                    b[k * ldb + iv + 1] = -b[(ki + 1) * ldb + iv + 1] * t[(ki + 1) * ldt + k]
                }
                var vmax = 1.0
                var vcrit = bignum
                var j = ki + 2
                while (j < n) {
                    if (j == n - 1 || t[(j + 1) * ldt + j] == 0.0) {
                        if (norms[j] > vcrit) {
                            let rec = 1.0 / vmax
                            nFloat64.dscal(n - ki, rec, b[ki * ldb + iv..], ldb)
                            nFloat64.dscal(n - ki, rec, b[ki * ldb + iv + 1..], ldb)
                            vmax = 1.0
                        }
                        b[j * ldb + iv] -= nFloat64.ddot(j - ki - 2, t[(ki + 2) * ldt + j..], ldt, b[(ki + 2) * ldb + iv..],
                            ldb)
                        b[j * ldb + iv + 1] -= nFloat64.ddot(j - ki - 2, t[(ki + 2) * ldt + j..], ldt,
                            b[(ki + 2) * ldb + iv + 1..], ldb)
                        let (scale, _, _) = dlaln2(false, 1, 2, smin, 1.0, t[j * ldt + j..], ldt, 1.0, 1.0,
                            b[j * ldb + iv..], ldb, wr, -wi, x[..2], 2)
                        if (scale != 1.0) {
                            nFloat64.dscal(n - ki, scale, b[ki * ldb + iv..], ldb)
                            nFloat64.dscal(n - ki, scale, b[ki * ldb + iv + 1..], ldb)
                        }
                        b[j * ldb + iv] = x[0]
                        b[j * ldb + iv + 1] = x[1]
                        vmax = nanMax(vmax, nanMax(abs(b[j * ldb + iv]), abs(b[j * ldb + iv + 1])))
                        vcrit = bignum / vmax
                        j++
                    } else {
                        if (nanMax(norms[j], norms[j + 1]) > vcrit) {
                            let rec = 1.0 / vmax
                            nFloat64.dscal(n - ki, rec, b[ki * ldb + iv..], ldb)
                            nFloat64.dscal(n - ki, rec, b[ki * ldb + iv + 1..], ldb)
                            vmax = 1.0
                        }
                        b[j * ldb + iv] -= nFloat64.ddot(j - ki - 2, t[(ki + 2) * ldt + j..], ldt, b[(ki + 2) * ldb + iv..],
                            ldb)
                        b[j * ldb + iv + 1] -= nFloat64.ddot(j - ki - 2, t[(ki + 2) * ldt + j..], ldt,
                            b[(ki + 2) * ldb + iv + 1..], ldb)
                        b[(j + 1) * ldb + iv] -= nFloat64.ddot(j - ki - 2, t[(ki + 2) * ldt + j + 1..], ldt,
                            b[(ki + 2) * ldb + iv..], ldb)
                        b[(j + 1) * ldb + iv + 1] -= nFloat64.ddot(j - ki - 2, t[(ki + 2) * ldt + j + 1..], ldt,
                            b[(ki + 2) * ldb + iv + 1..], ldb)
                        let (scale, _, _) = dlaln2(true, 2, 2, smin, 1.0, t[j * ldt + j..], ldt, 1.0, 1.0,
                            b[j * ldb + iv..], ldb, wr, -wi, x[..], 2)
                        if (scale != 1.0) {
                            nFloat64.dscal(n - ki, scale, b[ki * ldb + iv..], ldb)
                            nFloat64.dscal(n - ki, scale, b[ki * ldb + iv + 1..], ldb)
                        }
                        b[j * ldb + iv] = x[0]
                        b[j * ldb + iv + 1] = x[1]
                        b[(j + 1) * ldb + iv] = x[2]
                        b[(j + 1) * ldb + iv + 1] = x[3]
                        let vmax01 = nanMax(abs(x[0]), abs(x[1]))
                        let vmax23 = nanMax(abs(x[2]), abs(x[3]))
                        vmax = nanMax(vmax, nanMax(vmax01, vmax23))
                        vcrit = bignum / vmax
                        j += 2
                    }
                }
                if (howmny != EVAllMulQ) {
                    nFloat64.dcopy(n - ki, b[ki * ldb + iv..], ldb, vl[ki * ldvl + is_..], ldvl)
                    nFloat64.dcopy(n - ki, b[ki * ldb + iv + 1..], ldb, vl[ki * ldvl + is_ + 1..], ldvl)
                    var emax = 0.0
                    for (k in ki..n) {
                        emax = nanMax(emax, abs(vl[k * ldvl + is_]) + abs(vl[k * ldvl + is_ + 1]))
                    }
                    let remax = 1.0 / emax
                    nFloat64.dscal(n - ki, remax, vl[ki * ldvl + is_..], ldvl)
                    nFloat64.dscal(n - ki, remax, vl[ki * ldvl + is_ + 1..], ldvl)
                    for (k in 0..ki) {
                        vl[k * ldvl + is_] = 0.0
                        vl[k * ldvl + is_ + 1] = 0.0
                    }
                } else if (nb == 1) {
                    if (n - ki - 2 > 0) {
                        nFloat64.dgemv(NoTrans, n, n - ki - 2, 1.0, vl[ki + 2..], ldvl, b[(ki + 2) * ldb + iv..], ldb,
                            b[ki * ldb + iv], vl[ki..], ldvl)
                        nFloat64.dgemv(NoTrans, n, n - ki - 2, 1.0, vl[ki + 2..], ldvl, b[(ki + 2) * ldb + iv + 1..],
                            ldb, b[(ki + 1) * ldb + iv + 1], vl[ki + 1..], ldvl)
                    } else {
                        nFloat64.dscal(n, b[ki * ldb + iv], vl[ki..], ldvl)
                        nFloat64.dscal(n, b[(ki + 1) * ldb + iv + 1], vl[ki + 1..], ldvl)
                    }
                    var emax = 0.0
                    for (k in 0..n) {
                        emax = nanMax(emax, abs(vl[k * ldvl + ki]) + abs(vl[k * ldvl + ki + 1]))
                    }
                    let remax = 1.0 / emax
                    nFloat64.dscal(n, remax, vl[ki..], ldvl)
                    nFloat64.dscal(n, remax, vl[ki + 1..], ldvl)
                } else {
                    for (k in 0..ki) {
                        b[k * ldb + iv] = 0.0
                        b[k * ldb + iv + 1] = 0.0
                    }
                    iscomplex[iv] = ip
                    iscomplex[iv + 1] = -ip
                    iv++
                }
            }
            if (nb > 1) {
                var ki2 = ki
                if (ip != 0) {
                    ki2++
                }
                if (iv >= nb - 2 || ki2 == n - 1) {
                    nFloat64.dgemm(NoTrans, NoTrans, n, iv + 1, n - ki2 + iv, 1.0, vl[ki2 - iv..], ldvl,
                        b[(ki2 - iv) * ldb..], ldb, 0.0, b[nb..], ldb)
                    var remax: Float64 = 0.0
                    for (k in 0..=iv) {
                        if (iscomplex[k] == 0) {
                            let ii = nFloat64.idamax(n, b[nb + k..], ldb)
                            remax = 1.0 / abs(b[ii * ldb + nb + k])
                        } else if (iscomplex[k] == 1) {
                            var emax = 0.0
                            for (ii in 0..n) {
                                emax = nanMax(emax, abs(b[ii * ldb + nb + k]) + abs(b[ii * ldb + nb + k + 1]))
                            }
                            remax = 1.0 / emax
                        }
                        nFloat64.dscal(n, remax, b[nb + k..], ldb)
                    }
                    dlacpy(All, n, iv + 1, b[nb..], ldb, vl[ki2 - iv..], ldvl)
                    iv = 0
                } else {
                    iv++
                }
            }
            is_++
            if (ip != 0) {
                is_++
            }
        }

        return m
    }
}
