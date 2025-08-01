package cjnum.lapack.testlapack

import cjnum.blas.blas64.*
import cjnum.blas.*
import cjnum.floats.*
import cjnum.lapack.*
import cjnum.util.*
import std.math.*
import std.random.*
import std.sort.*

class dsterfTestStruct {
    dsterfTestStruct(
        let d!: Array<Float64>,
        let e!: Array<Float64>,
        let n!: Int64,
        let want!: Array<Float64>
    ) {}
}

public func dsterfTest(lFloat64: LFloat64): Bool {
    let tol = 1e-14

    var cas = 0
    for (test in [
        dsterfTestStruct(
            d: [1.0, 3.0, 4.0, 6.0],
            e: [2.0, 4.0, 5.0],
            n: 4,
            want: [11.046227528488854, 4.795922173417400, -2.546379458290125, 0.704229756383872]
        )
    ]) {
        let n = test.n
        let got = Array<Float64>(test.d.size, repeat: 0.0)
        copy(got, test.d)
        let e = Array<Float64>(test.e.size, repeat: 0.0)
        copy(e, test.e)
        let ok = lFloat64.dsterf(n, got, e)
        if (!ok) {
            throw Exception("Case ${cas}, n=${n}: Dsterf failed")
            // continue
        }
        let want = Array<Float64>(test.want.size, repeat: 0.0)
        copy(want, test.want)
        sort(want, stable: true)
        let diff = fDistance(Float64.Inf, s: got, t: want)
        if (diff > tol) {
            throw Exception("Case ${cas}, n=${n}: unexpected result, |dGot-dWant|=${diff}")
        }
        cas++
    }

    let rnd = Random(1)
    for (n in [0, 1, 2, 3, 4, 5, 6, 10, 50]) {
        for (typ in 0..=8) {
            let d = Array<Float64>(n, repeat: 0.0)
            let e = if (n > 1) {
                Array<Float64>(n - 1, repeat: 0.0)
            } else {
                Array<Float64>()
            }
            if (typ == 0) {
            } else if (typ == 1) {
                for (i in 0..d.size) {
                    d[i] = 1.0
                }
            } else if (typ == 2) {
                for (i in 0..n) {
                    if (i == 0) {
                        d[i] = 1.0
                    } else {
                        d[i] = 1.0 - (1.0 - dlamchE) * Float64(i) / Float64(n - 1)
                    }
                    if (rnd.nextFloat64() < 0.5) {
                        d[i] *= -1.0
                    }
                }
            } else if (typ == 3 || typ == 4 || typ == 5) {
                for (i in 0..n) {
                    if (i == 0) {
                        d[i] = 1.0
                    } else {
                        d[i] = pow(dlamchE, Float64(i) / Float64(n - 1))
                    }
                    if (rnd.nextFloat64() < 0.5) {
                        d[i] *= -1.0
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
                for (i in 0..d.size) {
                    if (rnd.nextFloat64() < 0.5) {
                        d[i] *= -1.0
                    }
                }
            } else if (typ == 7) {
                for (i in 0..d.size) {
                    d[i] = rnd.nextGaussianFloat64()
                }
            } else if (typ == 8) {
                for (i in 0..d.size) {
                    d[i] = rnd.nextGaussianFloat64()
                }
                for (i in 0..e.size) {
                    e[i] = rnd.nextGaussianFloat64()
                }
            }
            let eCopy = Array<Float64>(e.size, repeat: 0.0)
            copy(eCopy, e)

            let name = "n=${n},type=${typ}"

            let dGot = Array<Float64>(d.size, repeat: 0.0)
            copy(dGot, d)
            var ok = lFloat64.dsterf(n, dGot, e)
            if (!ok) {
                throw Exception("${name}: Dsterf failed")
                // continue
            }

            if (n == 0) {
                continue
            }

            var fis = true
            for (i in 0..dGot.size - 1) {
                if (dGot[i] > dGot[i + 1]) {
                    fis = false
                }
            }
            if (!fis) {
                throw Exception("${name}: eigenvalues are not sorted")
                // continue
            }

            let dWant = Array<Float64>(d.size, repeat: 0.0)
            copy(dWant, d)
            copy(e, eCopy)
            let z = nanGeneral(n, n, n)
            ok = lFloat64.dsteqr(EVTridiag, n, dWant, e, z.data, z.stride, Array<Float64>(2 * n, repeat: 0.0))
            if (!ok) {
                throw Exception("${name}: computing reference solution using Dsteqr failed")
                // continue
            }
            if (residualOrthogonal(z, false) > tol * Float64(n)) {
                throw Exception(
                    "${name}: Z is not orthogonal; resid=${residualOrthogonal(z, false)}, want<=${tol*Float64(n)}")
            }

            var diff = 0.0
            var dMax = 0.0
            for (i in 0..dGot.size) {
                let diffAbs = abs(dGot[i] - dWant[i])
                diff = nanMax(diff, diffAbs)
                let dAbs = nanMax(abs(dGot[i]), abs(dWant[i]))
                dMax = nanMax(dMax, dAbs)
            }
            dMax = nanMax(dlamchS, dMax)
            if (diff > tol * dMax) {
                throw Exception("${name}: unexpected result; |dGot-dWant|=${diff}")
            }

            copy(e, eCopy)
            let lda = n
            let a = Array<Float64>(n * lda, repeat: 0.0)
            var anorm = 0.0
            var tmp = 0.0
            for (i in 0..n - 1) {
                a[i * lda + i] = d[i]
                a[i * lda + i + 1] = e[i]
                let tmp2 = abs(e[i])
                anorm = nanMax(anorm, abs(d[i]) + tmp + tmp2)
                tmp = tmp2
            }
            a[(n - 1) * lda + n - 1] = d[n - 1]
            anorm = nanMax(anorm, abs(d[n - 1]) + tmp)

            let nFloat64 = nFloat64Implementation()
            for (i in 0..n) {
                nFloat64.dsyr(Upper, n, -dGot[i], z.data[i..], z.stride, a, lda)
            }

            let wnorm = lFloat64.dlansy(MaxColumnSum, Upper, n, a, lda, Array<Float64>(n, repeat: 0.0))

            if (anorm > wnorm) {
                diff = wnorm / anorm / Float64(n)
            } else {
                if (anorm < 1.0) {
                    diff = nanMin(wnorm, Float64(n) * anorm) / anorm / Float64(n)
                } else {
                    diff = nanMin(wnorm / anorm, Float64(n)) / Float64(n)
                }
            }

            if (diff > tol) {
                throw Exception("${name}: unexpected result; |A - Z D Zᵀ|/(|A| n)=${diff}")
            }
        }
    }
    return true
}
