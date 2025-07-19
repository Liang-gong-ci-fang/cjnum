package cjnum.lapack.cjnum

import cjnum.lapack.*
import std.math.*

extend LFloat64Implementation {
    public func dsterf(n: Int64, d: Array<Float64>, e: Array<Float64>): Bool {
        if (n < 0) {
            throw Exception(nLT0)
        }

        if (n == 0) {
            return true
        }

        if (d.size < n) {
            throw Exception(shortD)
        }
        if (e.size < n - 1) {
            throw Exception(shortE)
        }

        if (n == 1) {
            return true
        }

        let none = 0
        let down = 1
        let up = 2

        let eps = dlamchE
        let eps2 = eps * eps
        let safmin = dlamchS
        let safmax = 1.0 / safmin
        let ssfmax = sqrt(safmax) / 3.0
        let ssfmin = sqrt(safmin) / eps2

        let maxit = 30
        let nmaxit = n * maxit
        var jtot = 0
        var l1 = 0

        while (true) {
            if (l1 > n - 1) {
                dlasrt(SortIncreasing, n, d)
                return true
            }
            if (l1 > 0) {
                e[l1 - 1] = 0.0
            }
            var m = l1
            while (m < n - 1) {
                if (abs(e[m]) <= sqrt(abs(d[m])) * sqrt(abs(d[m + 1])) * eps) {
                    e[m] = 0.0
                    break
                }
                m++
            }

            var l = l1
            let lsv = l
            var lend = m
            let lendsv = lend
            l1 = m + 1
            if (lend == 0) {
                continue
            }

            let anorm = dlanst(MaxAbs, lend - l + 1, d[l..], e[l..])
            var iscale = none
            if (anorm == 0.0) {
                continue
            }
            if (anorm > ssfmax) {
                iscale = down
                dlascl(General, 0, 0, anorm, ssfmax, lend - l + 1, 1, d[l..], n)
                dlascl(General, 0, 0, anorm, ssfmax, lend - l, 1, e[l..], n)
            } else if (anorm < ssfmin) {
                iscale = up
                dlascl(General, 0, 0, anorm, ssfmin, lend - l + 1, 1, d[l..], n)
                dlascl(General, 0, 0, anorm, ssfmin, lend - l, 1, e[l..], n)
            }

            let el = e[l..lend]
            for (i in 0..el.size) {
                el[i] *= el[i]
            }

            if (abs(d[lend]) < abs(d[l])) {
                lend = lsv
                l = lendsv
            }
            if (lend >= l) {
                while (true) {
                    if (l != lend) {
                        m = l
                        while (m < lend) {
                            if (abs(e[m]) <= eps2 * abs(d[m] * d[m + 1])) {
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
                        (d[l], d[l + 1]) = dlae2(d[l], sqrt(e[l]), d[l + 1])
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

                    let rte = sqrt(e[l])
                    var sigma = (d[l+1] - p) / (2.0 * rte)
                    let r = dlapy2(sigma, 1.0)
                    sigma = p - (rte / (sigma + copysign(r, sigma)))

                    var c = 1.0
                    var s = 0.0
                    var gamma = d[m] - sigma
                    p = gamma * gamma

                    for (i in m - 1..=l: -1) {
                        let bb = e[i]
                        let r = p + bb
                        if (i != m - 1) {
                            e[i + 1] = s * r
                        }
                        let oldc = c
                        c = p / r
                        s = bb / r
                        let oldgam = gamma
                        let alpha = d[i]
                        gamma = c * (alpha - sigma) - s * oldgam
                        d[i + 1] = oldgam + alpha - gamma
                        if (c != 0.0) {
                            p = (gamma * gamma) / c
                        } else {
                            p = oldc * bb
                        }
                    }
                    e[l] = s * p
                    d[l] = sigma + gamma
                }
            } else {
                while (true) {
                    m = l
                    while (m > lend) {
                        if (abs(e[m - 1]) <= eps2 * abs(d[m] * d[m - 1])) {
                            break
                        }
                        m--
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
                        (d[l], d[l-1]) = dlae2(d[l], sqrt(e[l - 1]), d[l - 1])
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
                    let rte = sqrt(e[l - 1])
                    var sigma = (d[l-1] - p) / (2.0 * rte)
                    let r = dlapy2(sigma, 1.0)
                    sigma = p - (rte / (sigma + copysign(r, sigma)))

                    var c = 1.0
                    var s = 0.0
                    var gamma = d[m] - sigma
                    p = gamma * gamma

                    for (i in m..l) {
                        let bb = e[i]
                        let r = p + bb
                        if (i != m) {
                            e[i - 1] = s * r
                        }
                        let oldc = c
                        c = p / r
                        s = bb / r
                        let oldgam = gamma
                        let alpha = d[i + 1]
                        gamma = c * (alpha - sigma) - s * oldgam
                        d[i] = oldgam + alpha - gamma
                        if (c != 0.0) {
                            p = (gamma * gamma) / c
                        } else {
                            p = oldc * bb
                        }
                    }
                    e[l - 1] = s * p
                    d[l] = sigma + gamma
                }
            }

            if (iscale == down) {
                dlascl(General, 0, 0, ssfmax, anorm, lendsv - lsv + 1, 1, d[lsv..], n)
            }
            else if (iscale == up) {
                dlascl(General, 0, 0, ssfmin, anorm, lendsv - lsv + 1, 1, d[lsv..], n)
            }
            if (jtot >= nmaxit) {
                break
            }
        }
        for (v in e[0..n - 1]) {
            if (v != 0.0) {
                return false
            }
        }
        dlasrt(SortIncreasing, n, d)
        return true
    }
}