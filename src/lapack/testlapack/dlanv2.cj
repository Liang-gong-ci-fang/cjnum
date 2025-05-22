package cjnum.lapack.testlapack

import std.math.*
import std.random.*
import cjnum.lapack.*

public func dlanv2Test(lFloat64: LFloat64): Bool {
    let rnd = Random(1)
    for (_ in 0..10) {
        let a = rnd.nextFloat64()
        let b = rnd.nextFloat64()
        let d = rnd.nextFloat64()
        if (!_dlanv2Test(lFloat64, a, b, 0.0, d)) {
            return false
        }
    }
    for (i in 0..10) {
        let a = rnd.nextFloat64()
        let c = rnd.nextFloat64()
        let d = rnd.nextFloat64()
        if (!_dlanv2Test(lFloat64, a, 0.0, c, d)) {
            return false
        }
    }
    for (_ in 0..10) {
        let a = rnd.nextFloat64()
        let b = rnd.nextFloat64()
        var c = rnd.nextFloat64()
        if (b * c >= 0.0) {
            c = -c
        }
        if (!_dlanv2Test(lFloat64, a, b, c, a)) {
            return false
        }
    }
    for (i in 0..10) {
        let a = rnd.nextFloat64()
        let b = rnd.nextFloat64()
        let c = rnd.nextFloat64()
        let d = rnd.nextFloat64()
        if (!_dlanv2Test(lFloat64, a, b, c, d)) {
            return false
        }
    }
    _dlanv2Test(lFloat64, 0.0, 1.0, -1.0, Float64.MinNormal)
}

func _dlanv2Test(lFloat64: LFloat64, a: Float64, b: Float64, c: Float64, d: Float64): Bool {
    let aa: Float64
    let bb: Float64
    let cc: Float64
    let dd: Float64
    let cs: Float64
    let sn: Float64
    let rt1r: Float64
    let rt2r: Float64
    let rt1i: Float64
    let rt2i: Float64
    (aa, bb, cc, dd, rt1r, rt1i, rt2r, rt2i, cs, sn) = lFloat64.dlanv2(a, b, c, d)

    let mat = "[${a} ${b}; ${c} ${d}]"
    if (cc == 0.0) {
        if (rt1i != 0.0 || rt2i != 0.0) {
            throw Exception("Unexpected complex eigenvalues for ${mat}")
        }
    } else {
        if (aa != dd) {
            throw Exception("Diagonal elements not equal for ${mat}: got [${aa} ${dd}]")
        }
        if (bb * cc >= 0.0) {
            throw Exception("Non-diagonal elements have the same sign for ${mat}: got [${bb} ${cc}]")
        } else {
            let im = sqrt(-bb * cc)
            if (abs(rt1i - im) > 1e-14 && abs(rt1i + im) > 1e-14) {
                throw Exception("Unexpected imaginary part of eigenvalue for ${mat}: got ${rt1i}, want ${im} or ${-im}")
            }
            if (abs(rt2i - im) > 1e-14 && abs(rt2i + im) > 1e-14) {
                throw Exception("Unexpected imaginary part of eigenvalue for ${mat}: got ${rt2i}, want ${im} or ${-im}")
            }
        }
    }

    if (rt1r != aa && rt1r != dd) {
        throw Exception("Unexpected real part of eigenvalue for ${mat}: got ${rt1r}, want ${aa} or ${dd}")
    }
    if (rt2r != aa && rt2r != dd) {
        throw Exception("Unexpected real part of eigenvalue for ${mat}: got ${rt2r}, want ${aa} or ${dd}")
    }
    if (abs(hypot(cs, sn) - 1.0) > 1e-14) {
        throw Exception("Unexpected unitary matrix for ${mat}: got cs ${cs}, sn ${sn}")
    }

    let gota = cs * (aa * cs - bb * sn) - sn * (cc * cs - dd * sn)
    let gotb = cs * (aa * sn + bb * cs) - sn * (cc * sn + dd * cs)
    let gotc = sn * (aa * cs - bb * sn) + cs * (cc * cs - dd * sn)
    let gotd = sn * (aa * sn + bb * cs) + cs * (cc * sn + dd * cs)
    if (abs(gota - a) > 1e-14 || abs(gotb - b) > 1e-14 || abs(gotc - c) > 1e-14 || abs(gotd - d) > 1e-14) {
        throw Exception(
            "Unexpected factorization: got [${gota}, ${gotb}, ${gotc}, ${gotd}], want [${a}, ${b}, ${c}, ${d}]")
    }
    return true
}
