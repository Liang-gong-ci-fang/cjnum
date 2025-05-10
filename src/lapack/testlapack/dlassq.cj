package cjnum.lapack.testlapack

import std.math.*
import std.random.*
import cjnum.floats.*
import cjnum.lapack.*

class dlassqTestStruct {
    dlassqTestStruct(
        let x!: Array<Float64>,
        let a!: Float64
    ) { }
}

public func dlassqTest(lFloat64: LFloat64): Bool {
    let values = [
		0.0,
		0.5 * safmin,
		smlnum,
		ulp,
		1.0,
		1.0 / ulp,
		bignum,
		safmax,
		Float64.Inf,
		Float64.NaN
    ]

    let rnd = Random(1)
    for (n in [0, 1, 2, 3, 4, 5, 10, 20, 30, 40]) {
        for (incx in [1, 3]) {
            for (cas in 0..3) {
                for (v0_ in values) {
                    var v0 = v0_
                    if (v0 > 1.0) {
                        v0 *= 0.5
                    }
                    for (v1_ in values) {
                        var v1 = v1_
                        if (v1 > 1.0) {
                            v1 = 0.5 * v1 / sqrt(Float64(n + 1))
                        }
                        _dlassqTest(lFloat64, rnd, n, incx, cas, v0, v1)
                    }
                }
            }
        }
    }
    return true
}

func _dlassqTest(lFloat64: LFloat64, rnd: Random, n: Int64, incx: Int64, cas: Int64, v0: Float64, v1: Float64): Unit {
    const rogue = 1234.5678
    const tol = 1e-15

    let name = "n = ${n}, incx = ${incx}, cas = ${cas}, v0 = ${v0}, v1 = ${v1}"

    let work = Array<Float64>(n, { _ => 
            1.0 - 2.0 * rnd.nextFloat64()
        }
    )

    var workssq = 0.0
    for (wi in work) {
        workssq += wi * wi
    }

    var scale = 0.0
    var ssq = 0.0
    if (cas == 0) {
        scale = 1.0
        ssq = v0
    } else if (cas == 1) {
        scale = sqrt(v0)
        ssq = 1.0
    } else if (cas == 2) {
        if (v0 < 1.0) {
            scale = 1.0 / 3.0
            ssq = 9.0 * v0
        } else {
            scale = 3.0
            ssq = v0 / 9.0
        }
    } else {
        throw Exception("bad cas")
    }

    let x = Array<Float64>(max(0, 1 + (n - 1) * incx), item: rogue)
    for (i in 0..work.size) {
        x[i * incx] = v1 * work[i]
    }
    let xCopy = x.clone()

    let (scaleGot, ssqGot) = lFloat64.dlassq(n, x, incx, scale, ssq)
    let nrmGot = scaleGot * sqrt(ssqGot)

    if (!fSame(x, xCopy)) {
        throw Exception("${name}: unexpected modification of x")
    }

    let z0 = sqrt(v0)
    var z1n = 0.0
    if (n >= 1) {
        z1n = v1 * sqrt(workssq)
    }
    let zmin = min(z0, z1n)
    let zmax = max(z0, z1n)
    var nrmWant = 0.0

    if (z0.isNaN() || z1n.isNaN()) {
        nrmWant = Float64.NaN
    } else if (zmin == zmax) {
        nrmWant = sqrt(2.0) * zmax
    } else if (zmax == 0.0) {
        nrmWant = 0.0
    } else {
        nrmWant = zmax * sqrt(1.0 + (zmin / zmax) * (zmin / zmax))
    }

    if (nrmGot.isNaN() || nrmWant.isNaN()) {
        if (!nrmGot.isNaN()) {
            throw Exception("${name}: expected NaN; got ${nrmGot}")
        }
        if (!nrmWant.isNaN()) {
            throw Exception("${name}: unexpected NaN; want ${nrmWant}")
        }
    } else if (nrmGot == nrmWant) {
        // do nothing
    } else if (nrmWant == 0.0) {
        if (nrmGot > tol) {
            throw Exception("${name}: unexpected result; got ${nrmGot}, want 0")
        }
    } else {
        let diff = abs(nrmGot - nrmWant) / nrmWant / max(1.0, Float64(n))
        if (diff.isNaN() || diff > tol) {
            throw Exception("${name}: unexpected result; got ${nrmGot}, want ${nrmWant}")
        }
    }
}