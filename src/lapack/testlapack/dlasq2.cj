package cjnum.lapack.testlapack

import cjnum.blas.*
import cjnum.floats.*
import cjnum.lapack.*
import std.math.*
import std.random.*
import std.sort.*

public func dlasq2Test(lFloat64: LFloat64): Bool {
    let tol = 1e-14
    let rnd = Random(1)

    for (n in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 25, 50]) {
        for (_ in 0..10) {
            for (typ in 0..=2) {
                let name = "n=${n},typ=${typ}"

                let want = Array<Float64>(n, item: 0.0)
                let z = Array<Float64>(4 * n, item: 0.0)

                if (typ == 0) {
                } else if (typ == 1) {
                    for (i in 0..n) {
                        z[2 * i] = rnd.nextFloat64()
                        want[i] = z[2 * i]
                    }
                    stableSort(want)
                } else if (typ == 2) {
                    for (i in 0..z.size) {
                        z[i] = rnd.nextFloat64()
                    }
                }

                let zCopy = Array<Float64>(z.size, item: 0.0)
                copy(zCopy, z)

                let info = lFloat64.dlasq2(n, z)
                if (info != 0) {
                    throw Exception("${name}: Dlasq2 failed")
                }
                if (n == 0) {
                    continue
                }

                let got = z[..n]

                if (typ == 2) {
                    let ldt = n
                    let T = Array<Float64>(n * ldt, item: 0.0)

                    for (i in 0..n) {
                        if (i == 0) {
                            T[0] = zCopy[0]
                        } else {
                            T[i * ldt + i] = zCopy[2 * i - 1] + zCopy[2 * i]
                        }
                        if (i < n - 1) {
                            T[i * ldt + i + 1] = 1.0
                            T[(i + 1) * ldt + i] = zCopy[2 * i + 1] * zCopy[2 * i]
                        }
                    }

                    let d = Array<Float64>(n, item: 0.0)
                    d[0] = 1.0
                    for (i in 1..n) {
                        d[i] = d[i - 1] * T[i * ldt + i - 1] / T[(i - 1) * ldt + i]
                    }
                    for (i in 0..n) {
                        d[i] = sqrt(d[i])
                    }
                    for (i in 0..n) {
                        for (j in i..=min(i + 1, n - 1)) {
                            T[i * ldt + j] *= d[j] / d[i]
                        }
                    }

                    let work = Array<Float64>(3 * n, item: 0.0)
                    let ok = lFloat64.dsyev(EVNone, Upper, n, T, ldt, want, work, work.size)
                    if (!ok) {
                        throw Exception("${name}: Dsyev failed")
                    }
                }

                stableSort(got)
                let diff = fDistance(Float64.Inf, s: got, t: want)
                if (diff > tol) {
                    throw Exception("${name}: unexpected eigenvalues; diff=${diff}\n${got}\n%{want}\n\n")
                }
            }
        }
    }
    return true
}
