package cjnum.lapack.cjnum

import cjnum.blas.blas64.*
import cjnum.blas.*
import cjnum.lapack.General as lGeneral
import cjnum.lapack.*
import std.math.*

extend LFloat64Implementation {
    public func dsteqr(compz: EVComp, n: Int64, d: Array<Float64>, e: Array<Float64>, z: Array<Float64>, ldz: Int64,
        work: Array<Float64>): Bool {
        if (compz != EVCompNone && compz != EVTridiag && compz != EVOrig) {
            throw Exception(badEVComp)
        } else if (n < 0) {
            throw Exception(nLT0)
        } else if (ldz < 1 || (compz != EVCompNone && ldz < n)) {
            throw Exception(badLdZ)
        }

        if (n == 0) {
            return true
        }

        if (d.size < n) {
            throw Exception(shortD)
        } else if (e.size < n - 1) {
            throw Exception(shortE)
        } else if (compz != EVCompNone && z.size < (n - 1) * ldz + n) {
            throw Exception(shortZ)
        } else if (compz != EVCompNone && work.size < max(1, 2 * n - 2)) {
            throw Exception(shortWork)
        }

        var icompz = 0
        if (compz == EVOrig) {
            icompz = 1
        } else if (compz == EVTridiag) {
            icompz = 2
        }

        if (n == 1) {
            if (icompz == 2) {
                z[0] = 1.0
            }
            return true
        }

        let nFloat64 = nFloat64Implementation()

        let eps = dlamchE
        let eps2 = eps * eps
        let safmin = dlamchS
        let safmax = 1.0 / safmin
        let ssfmax = sqrt(safmax) / 3.0
        let ssfmin = sqrt(safmin) / eps2

        if (icompz == 2) {
            dlaset(All, n, n, 0.0, 1.0, z, ldz)
        }
        const maxit = 30
        let nmaxit = n * maxit

        var jtot = 0

        var l1 = 0
        let nm1 = n - 1

        const down = 1
        const up = 2
        var iscale = 0

        while (true) {
            if (l1 > n - 1) {
                if (icompz == 0) {
                    dlasrt(SortIncreasing, n, d)
                } else {
                    for (ii in 1..n) {
                        let i = ii - 1
                        var k = i
                        var p = d[i]
                        for (j in ii..n) {
                            if (d[j] < p) {
                                k = j
                                p = d[j]
                            }
                        }
                        if (k != i) {
                            d[k] = d[i]
                            d[i] = p
                            nFloat64.dswap(n, z[i..], ldz, z[k..], ldz)
                        }
                    }
                }
                return true
            }
            if (l1 > 0) {
                e[l1 - 1] = 0.0
            }
            var m = 0
            if (l1 <= nm1) {
                m = l1
                while (m < nm1) {
                    let test = abs(e[m])
                    if (test == 0.0) {
                        break
                    }
                    if (test <= (sqrt(abs(d[m])) * sqrt(abs(d[m + 1]))) * eps) {
                        e[m] = 0.0
                        break
                    }
                    m++
                }
            }
            var l = l1
            let lsv = l
            var lend = m
            let lendsv = lend
            l1 = m + 1
            if (lend == l) {
                continue
            }

            let anorm = dlanst(MaxAbs, lend - l + 1, d[l..], e[l..])
            if (anorm == 0.0) {
                continue
            } else if (anorm > ssfmax) {
                iscale = down
                dlascl(lGeneral, 0, 0, anorm, ssfmax, lend - l + 1, 1, d[l..], 1)
                dlascl(lGeneral, 0, 0, anorm, ssfmax, lend - l, 1, e[l..], 1)
            } else if (anorm < ssfmin) {
                iscale = up
                dlascl(lGeneral, 0, 0, anorm, ssfmin, lend - l + 1, 1, d[l..], 1)
                dlascl(lGeneral, 0, 0, anorm, ssfmin, lend - l, 1, e[l..], 1)
            }

            if (abs(d[lend]) < abs(d[l])) {
                lend = lsv
                l = lendsv
            }
            if (lend > l) {
                while (true) {
                    if (l != lend) {
                        m = l
                        while (m < lend) {
                            let v = abs(e[m])
                            if (v * v <= (eps2 * abs(d[m])) * abs(d[m + 1]) + safmin) {
                                break
                            }
                            m++
                        }
                    } else {
                        m = lend
                    }
                    if (m < lend) {
                        e[m] = 0.0
                    }
                    var p = d[l]
                    if (m == l) {
                        l++
                        if (l > lend) {
                            break
                        }
                        continue
                    }

                    if (m == l + 1) {
                        if (icompz > 0) {
                            (d[l], d[l + 1], work[l], work[n - 1 + l]) = dlaev2(d[l], e[l], d[l + 1])
                            dlasr(Right, Variable, Backward, n, 2, work[l..], work[n - 1 + l..], z[l..], ldz)
                        } else {
                            (d[l], d[l + 1]) = dlae2(d[l], e[l], d[l + 1])
                        }
                        e[l] = 0.0
                        l += 2
                        if (l > lend) {
                            break
                        }
                        continue
                    }

                    if (jtot == nmaxit) {
                        break
                    }
                    jtot++

                    var g = (d[l + 1] - p) / (2.0 * e[l])
                    var r = dlapy2(g, 1.0)
                    g = d[m] - p + e[l] / (g + copysign(r, g))
                    var s = 1.0
                    var c = 1.0
                    p = 0.0

                    for (i in m - 1..=l : -1) {
                        let f = s * e[i]
                        let b = c * e[i]
                        (c, s, r) = dlartg(g, f)
                        if (i != m - 1) {
                            e[i + 1] = r
                        }
                        g = d[i + 1] - p
                        r = (d[i] - g) * s + 2.0 * c * b
                        p = s * r
                        d[i + 1] = g + p
                        g = c * r - b

                        if (icompz > 0) {
                            work[i] = c
                            work[n - 1 + i] = -s
                        }
                    }
                    if (icompz > 0) {
                        let mm = m - l + 1
                        dlasr(Right, Variable, Backward, n, mm, work[l..], work[n - 1 + l..], z[l..], ldz)
                    }
                    d[l] -= p
                    e[l] = g
                }
            } else {
                while (true) {
                    if (l != lend) {
                        m = l
                        while (m > lend) {
                            let v = abs(e[m - 1])
                            if (v * v <= (eps2 * abs(d[m]) * abs(d[m - 1]) + safmin)) {
                                break
                            }
                            m--
                        }
                    } else {
                        m = lend
                    }
                    if (m > lend) {
                        e[m - 1] = 0.0
                    }
                    var p = d[l]
                    if (m == l) {
                        l--
                        if (l < lend) {
                            break
                        }
                        continue
                    }

                    if (m == l - 1) {
                        if (icompz > 0) {
                            (d[l - 1], d[l], work[m], work[n - 1 + m]) = dlaev2(d[l - 1], e[l - 1], d[l])
                            dlasr(Right, Variable, Forward, n, 2, work[m..], work[n - 1 + m..], z[l - 1..], ldz)
                        } else {
                            (d[l - 1], d[l]) = dlae2(d[l - 1], e[l - 1], d[l])
                        }
                        e[l - 1] = 0.0
                        l -= 2
                        if (l < lend) {
                            break
                        }
                        continue
                    }
                    if (jtot == nmaxit) {
                        break
                    }
                    jtot++

                    var g = (d[l - 1] - p) / (2.0 * e[l - 1])
                    var r = dlapy2(g, 1.0)
                    g = d[m] - p + (e[l - 1]) / (g + copysign(r, g))
                    var s = 1.0
                    var c = 1.0
                    p = 0.0

                    for (i in m..l) {
                        let f = s * e[i]
                        let b = c * e[i]
                        (c, s, r) = dlartg(g, f)
                        if (i != m) {
                            e[i - 1] = r
                        }
                        g = d[i] - p
                        r = (d[i + 1] - g) * s + 2.0 * c * b
                        p = s * r
                        d[i] = g + p
                        g = c * r - b

                        if (icompz > 0) {
                            work[i] = c
                            work[n - 1 + i] = s
                        }
                    }

                    if (icompz > 0) {
                        let mm = l - m + 1
                        dlasr(Right, Variable, Forward, n, mm, work[m..], work[n - 1 + m..], z[m..], ldz)
                    }
                    d[l] -= p
                    e[l - 1] = g
                }
            }

            if (iscale == down) {
                dlascl(lGeneral, 0, 0, ssfmax, anorm, lendsv - lsv + 1, 1, d[lsv..], 1)
                dlascl(lGeneral, 0, 0, ssfmax, anorm, lendsv - lsv, 1, e[lsv..], 1)
            } else if (iscale == up) {
                dlascl(lGeneral, 0, 0, ssfmin, anorm, lendsv - lsv + 1, 1, d[lsv..], 1)
                dlascl(lGeneral, 0, 0, ssfmin, anorm, lendsv - lsv, 1, e[lsv..], 1)
            }

            if (jtot >= nmaxit) {
                break
            }
        }
        for (i in 0..n - 1) {
            if (e[i] != 0.0) {
                return false
            }
        }
        return true
    }
}
