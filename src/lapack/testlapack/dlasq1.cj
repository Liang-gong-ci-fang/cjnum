package cjnum.lapack.testlapack

import cjnum.floats.*
import cjnum.lapack.*
import std.collection.*
import std.math.*
import std.random.*
import std.sort.*

public func dlasq1Test(lFloat64: LFloat64): Bool {
    const tol = 1e-14
    let rnd = Random(1)
    for (n in [0, 1, 2, 3, 4, 5, 8, 10, 30, 50]) {
        for (typ in 0..=7) {
            let name = "n=${n},typ=${typ}"

            let d = Array<Float64>(n, repeat: 0.0)
            if (typ == 0) {
            } else if (typ == 1) {
                for (i in 0..d.size) {
                    d[i] = 1.0
                }
            } else if (typ == 2) {
                for (i in 0..n) {
                    if (i == 0) {
                        d[0] = 1.0
                    } else {
                        d[i] = 1.0 - (1.0 - dlamchE) * Float64(i) / Float64(n - 1)
                    }
                }
            } else if (typ == 3 || typ == 4 || typ == 5) {
                for (i in 0..n) {
                    if (i == 0) {
                        d[0] = 1.0
                    } else {
                        d[i] = pow(dlamchE, Float64(i) / Float64(n - 1))
                    }
                }
                if (typ == 4) {
                    fScale(sqrt(1.0 / dlamchS), d)
                } else if (typ == 5) {
                    fScale(sqrt(dlamchS), d)
                }
            } else if (typ == 6) {
                for (i in 0..d.size) {
                    if (i == 0) {
                        d[i] = 1.0
                    } else {
                        d[i] = dlamchE
                    }
                }
            } else if (typ == 7) {
                for (i in 0..d.size) {
                    d[i] = abs(rnd.nextGaussianFloat64())
                }
            }

            let dWant = Array<Float64>(n, repeat: 0.0)
            copy(dWant, d)
            // sort.Sort(sort.Reverse(sort.Float64Slice(dWant)))
            sort(dWant)
            dWant.reverse()

            let work = Array<Float64>(max(1, 4 * n), repeat: 0.0)

            let lda = max(1, n)
            let a = Array<Float64>(n * lda, repeat: 0.0)
            dlagge(n, n, 0, 0, d, a, lda, rnd, work)

            let tauQ = Array<Float64>(n, repeat: 0.0)
            let tauP = Array<Float64>(n, repeat: 0.0)
            let e = Array<Float64>(max(0, n - 1), repeat: 0.0)
            lFloat64.dgebrd(n, n, a, lda, d, e, tauQ, tauP, work, work.size)

            for (i in 0..work.size) {
                work[i] = Float64.NaN
            }
            let info = lFloat64.dlasq1(n, d, e, work)
            if (info != 0) {
                // t.Fatalf("${name}: Dlasq1 returned non-zero info=${info}")
                throw Exception("${name}: Dlasq1 returned non-zero info=${info}")
            }

            if (n == 0) {
                continue
            }
            
            var srf = true
            for (i in 0..d.size - 1) {
                if (d[i] > d[i + 1]) {
                    srf = false
                }
            }

            if (!srf) {
                // t.Errorf("${name}: singular values not sorted")
                throw Exception("${name}: singular values not sorted")
            }

            let diff = fDistance(Float64.Inf, s: d, t: dWant)
            let md = match (max(dWant)) {
                case Some(v) => v
                case None => 0.0
            }
            if (diff > tol * md) {
                // t.Errorf("${name}: unexpected result; diff=${diff}")
                throw Exception("${name}: unexpected result; diff=${diff}")
            }
        }
    }
    return true
}
