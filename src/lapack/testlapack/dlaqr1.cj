package cjnum.lapack.testlapack

import std.math.*
import std.random.*
import cjnum.floats.*
import cjnum.lapack.*
import cjnum.blas.*
import cjnum.complex.*

public func dlaqr1Test(lFloat64: LFloat64): Bool {
    // let rnd = Random(1)

    for (n in [2, 3]) {
        let rnd = Random(1)
        for (ldh in [n, n + 1, n + 10]) {
            for (cas in [1, 2]) {
                for (k in 0..100) {
                    let v = Array<Float64>(n, { _ => Float64.NaN })
                    let h = Array<Float64>(n * (n - 1) * ldh, { _ => Float64.NaN })
                    for (i in 0..n) {
                        for (j in 0..n) {
                            h[i * ldh + j] = rnd.nextGaussianFloat64()
                        }
                    }
                    var sr1 = 0.0
                    var sr2 = 0.0
                    var si1 = 0.0
                    var si2 = 0.0
                    if (cas == 1) {
                        sr1 = rnd.nextGaussianFloat64()
                        sr2 = sr1
                        si1 = rnd.nextGaussianFloat64()
                        si2 = -si1
                    } else {
                        sr1 = rnd.nextGaussianFloat64()
                        sr2 = rnd.nextGaussianFloat64()
                        si1 = 0.0
                        si2 = 0.0
                    }

                    lFloat64.dlaqr1(n, h, ldh, sr1, si1, sr2, si2, v)

                    let h1 = Array<Complex128>(n * n, { _ => Complex128(0.0, 0.0) })
                    for (i in 0..n) {
                        for (j in 0..n) {
                            h1[i * n + j] = Complex128(h[i * ldh + j], 0.0)
                            if (i == j) {
                                h1[i * n + j] -= Complex128(sr1, si1)
                            }
                        }
                    }

                    let h2 = Array<Complex128>(n, { i => Complex128(h[i * ldh], 0.0) })
                    h2[0] -= Complex128(sr2, si2)
                    let wantv = Array<Float64>(n, repeat: 0.0)
                    for (i in 0..n) {
                        for (j in 0..n) {
                            wantv[i] += real(h1[i * n + j] * h2[j])
                        }
                    }

                    let scale = v[0] / wantv[0]
                    fScale(scale, wantv)
                    if (fDistance(Float64.Inf, s: wantv, t: v) > 1e-13) {
                        throw Exception("n = ${n}, ldh = ${ldh}, case = ${cas}, k = ${k}: Unexpected value of v: got ${v}, want ${wantv}")
                    }
                }
            }
        }
    }
    return true
}