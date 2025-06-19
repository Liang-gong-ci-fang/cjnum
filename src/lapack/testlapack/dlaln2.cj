package cjnum.lapack.testlapack

import cjnum.complex.*
import cjnum.lapack.*
import std.math.*
import std.random.*

public func dlaln2Test(lFloat64: LFloat64): Bool {
    let rnd = Random(1)
    for (trans in [true, false]) {
        for (na in [1, 2]) {
            for (nw in [1, 2]) {
                for (extra in [0, 1, 2, 13]) {
                    for (_ in 0..1000) {
                        if (!testDlaln2(lFloat64, trans, na, nw, extra, rnd)) {
                            return false
                        }
                    }
                }
            }
        }
    }
    return true
}

func testDlaln2(lFloat64: LFloat64, trans: Bool, na: Int64, nw: Int64, extra: Int64, rnd: Random): Bool {
    const tol = 1e-11

    let ca = rnd.nextGaussianFloat64()
    let d1 = rnd.nextGaussianFloat64()
    let d2 = rnd.nextGaussianFloat64()

    var w: Complex128
    if (nw == 1) {
        w = Complex128(rnd.nextGaussianFloat64(), 0.0)
    } else {
        w = Complex128(rnd.nextGaussianFloat64(), rnd.nextGaussianFloat64())
    }
    let smin = dlamchP * (abs(w.real) + abs(w.imag))

    let a = randomGeneral(na, na, na + extra, rnd)
    let b = randomGeneral(na, nw, nw + extra, rnd)
    let x = randomGeneral(na, nw, nw + extra, rnd)

    let (scale, xnormGot, ok) = lFloat64.dlaln2(trans, na, nw, smin, ca, a.data, a.stride, d1, d2, b.data, b.stride,
        w.real, w.imag, x.data, x.stride)

    let prefix = "n=${na},lda=${a.stride},big=${b.stride}"

    if (!generalOutsideAllNaN(a)) {
        throw Exception("${prefix}: out-of-range write to A\n${a.data}")
    }
    if (!generalOutsideAllNaN(b)) {
        throw Exception("${prefix}: out-of-range write to B\n${b.data}")
    }
    if (!generalOutsideAllNaN(x)) {
        throw Exception("${prefix}: out-of-range write to X\n${x.data}")
    }

    if (scale <= 0.0 || 1.0 < scale) {
        throw Exception("${prefix}: invalid value of scale=${scale}")
    }

    var xnormWant: Float64 = 0.0
    for (i in 0..na) {
        var rowsum: Float64 = 0.0
        for (j in 0..nw) {
            rowsum += abs(x.data[i * x.stride + j])
        }
        if (rowsum > xnormWant) {
            xnormWant = rowsum
        }
    }
    if (xnormWant != xnormGot) {
        throw Exception("Case ${prefix}: unexpected xnorm with scale=${scale}. Want ${xnormWant}, got ${xnormGot}")
    }

    if (!ok) {
        return true
    }

    let m = Array<Complex128>(na * na, item: Complex128(0.0))
    if (trans) {
        for (i in 0..na) {
            for (j in 0..na) {
                m[i * na + j] = Complex128(ca * a.data[j * a.stride + i], 0.0)
            }
        }
    } else {
        for (i in 0..na) {
            for (j in 0..na) {
                m[i * na + j] = Complex128(ca * a.data[i * a.stride + j], 0.0)
            }
        }
    }
    m[0] -= w * Complex128(d1, 0.0)
    if (na == 2) {
        m[3] -= w * Complex128(d2, 0.0)
    }

    let cx = Array<Complex128>(na, item: Complex128(0.0))
    let cb = Array<Complex128>(na, item: Complex128(0.0))

    if (nw == 1) {
        for (i in 0..na) {
            cx[i] = Complex128(x.data[i * x.stride], 0.0)
            cb[i] = Complex128(scale * b.data[i * x.stride], 0.0)
        }
    } else if (nw == 2) {
        for (i in 0..na) {
            cx[i] = Complex128(x.data[i * x.stride], x.data[i * x.stride + 1])
            cb[i] = Complex128(scale * b.data[i * b.stride], scale * b.data[i * b.stride + 1])
        }
    }

    let mx = Array<Complex128>(na, item: Complex128(0.0))
    for (i in 0..na) {
        for (j in 0..na) {
            mx[i] += m[i * na + j] * cx[j]
        }
    }
    for (i in 0..na) {
        if (c128Abs(mx[i] - cb[i]) > tol) {
            throw Exception(
                "Case ${prefix}: unexpected value of left-hand side at row ${i} with scale=${scale}. Want ${cb[i]}, got ${mx[i]}"
            )
        }
    }
    return true
}
