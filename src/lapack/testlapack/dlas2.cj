package cjnum.lapack.testlapack

import std.math.*
import cjnum.lapack.*

class dlas2TestStruct {
    dlas2TestStruct(
        let f: Float64,
        let g: Float64,
        let h: Float64,
        let ssmin: Float64,
        let ssmax: Float64
    ) {}
}

public func dlas2Test(lFloat64: LFloat64): Bool {
    let tests = [
		dlas2TestStruct(10.0, 30.0, 12.0, 3.567778859365365, 33.634371616111189),
		dlas2TestStruct(10.0, 30.0, -12.0, 3.567778859365365, 33.634371616111189),
		dlas2TestStruct(2.0, 30.0, -12.0, 0.741557056404952, 32.364333658088754),
		dlas2TestStruct(-2.0, 5.0, 12.0, 1.842864429909778, 13.023204317408728)
    ]
    for (i in 0..tests.size) {
        let test = tests[i]
        let (ssmin, ssmax) = lFloat64.dlas2(test.f, test.g, test.h)
        if (abs(ssmin - test.ssmin) > 1e-12) {
            throw Exception("Case ${i}, minimal singular value mismatch. Want ${test.ssmin}, got ${ssmin}")
        }
        if (abs(ssmax - test.ssmax) > 1e-12) {
            throw Exception("Case ${i}, minimal singular value mismatch. Want ${test.ssmax}, got ${ssmax}")
        }
    }
    return true
}