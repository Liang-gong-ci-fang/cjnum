package cjnum.lapack.testlapack

import std.math.*
import std.random.*
import cjnum.lapack.*
import cjnum.util.*

public func dpttrfTest(lFloat64: LFloat64): Bool {
    let rnd = Random(1)
    for (n in [0, 1, 2, 3, 4, 5, 10, 20, 50, 51, 52, 53, 54, 100]) {
        _dpttrfTest(lFloat64, rnd, n)
    }
    return true
}

func _dpttrfTest(lFloat64: LFloat64, rnd: Random, n: Int64): Unit {
	const tol = 1e-15

	let name = "n=${n}"

	let (d, e) = newRandomSymTridiag(n, rnd)

	var dFac: Array<Float64> = [] 
    var eFac: Array<Float64> = []
	if (n > 0) {
		dFac = Array<Float64>(d.size, repeat: 0.0)
		copy(dFac, d)
		if (n > 1) {
			eFac = Array<Float64>(e.size, repeat: 0.0)
			copy(eFac, e)
		}
	}

	let ok = lFloat64.dpttrf(n, dFac, eFac)
	if (!ok) {
		throw Exception("${name}: bad test matrix, Dpttrf failed")
	}

	let resid = dpttrfResidual(n, d, e, dFac, eFac)
	if (resid > tol) {
		throw Exception("${name}: unexpected residual |L*D*Lᵀ - A|/(n * norm(A)); got ${resid}, want <= ${tol}")
	}
}

func dpttrfResidual(n: Int64, d: Array<Float64>, e: Array<Float64>, dFac: Array<Float64>, eFac: Array<Float64>): Float64 {
    if (n == 0) {
        return 0.0
    }

    let dDiff = Array<Float64>(n, repeat: 0.0)
    let eDiff = Array<Float64>(n - 1, repeat: 0.0)
    dDiff[0] = dFac[0] - d[0]
    for (i in 0..eFac.size) {
        let ef = eFac[i]
        let de = dFac[i] * ef
		dDiff[i + 1] = de * ef + dFac[i + 1] - d[i + 1]
		eDiff[i] = de - e[i]
    }

    var resid = 0.0
    if (n == 1) {
        resid = abs(dDiff[0])
    } else {
        resid = nanMax(abs(dDiff[0]) + abs(eDiff[0]), abs(dDiff[n - 1]) + abs(eDiff[n - 2]))
        for (i in 1..n - 1) {
            resid = nanMax(resid, abs(dDiff[i]) + abs(eDiff[i - 1]) + abs(eDiff[i]))
        }
    }

    let anorm = dlanst(MaxColumnSum, n, d, e)
    if (anorm == 0.0) {
        if (resid != 0.0) {
            return Float64.Inf
        }
        return 0.0
    }
    return resid / Float64(n) / anorm
}