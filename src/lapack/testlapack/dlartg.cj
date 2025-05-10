package cjnum.lapack.testlapack

import std.math.*
import cjnum.lapack.*

public func dlartgTest(lFloat64: LFloat64): Bool {
    const tol = 20.0 * ulp

    let values = [
		-safmax,
		-1.0 / ulp,
		-1.0,
		-1.0 / 3.0,
		-ulp,
		-safmin,
		0.0,
		safmin,
		ulp,
		1.0 / 3.0,
		1.0,
		1.0 / ulp,
		safmax,
		-Float64.Inf,
		Float64.Inf,
		Float64.NaN
    ]

    for (f in values) {
        for (g in values) {
            let name = "Case f=${f},g=${g}"

            let (cs, sn, r) = lFloat64.dlartg(f, g)

            if (f.isNaN() || g.isNaN()) {
                if (!r.isNaN()) {
                    throw Exception("${name}: unexpected r=${r}; want NaN")
                }
            } else if (f.isInf() || g.isInf()) {
                if (!r.isNaN() && !r.isInf()) {
                    throw Exception("${name}: unexpected r=${r}; want NaN or Inf")
                }
            } else {
                var d = max(abs(f), abs(g))
                d = min(max(safmin, d), safmax)
                let fs = f / d
                let gs = g / d
                let rs = r / d

                var rnorm = abs(rs)
                if (rnorm == 0.0) {
                    rnorm = max(abs(fs), abs(gs))
                    if (rnorm == 0.0) {
                        rnorm == 1.0
                    }
                }
                var resid = abs(rs - (cs * fs + sn * gs)) / rnorm
                if (resid > tol) {
                    throw Exception("${name}: cs*f + sn*g != r; resid=${resid}")
                }

                resid = abs(-sn * fs + cs * gs)
                if (resid > tol) {
                    throw Exception("${name}: -sn*f + cs*g != 0; resid=${resid}")
                }

                resid = abs(1.0 - (cs * cs + sn * sn))
                if (resid > tol) {
                    throw Exception("${name}: cs*cs + sn*sn != 1; resid=${resid}")
                }

                if (abs(f) > abs(g) && cs < 0.0) {
                    throw Exception("${name}: cs is negative; cs=${cs}")
                }
            }
        }
    }
    return true
}
