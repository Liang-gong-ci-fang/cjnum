package cjnum.lapack.testlapack

import std.math.*
import cjnum.lapack.*

public func dlae2Test(lFloat64: LFloat64): Bool {
    for (test in [
		[-10.0, 5.0, 3.0],
		[3.0, 5.0, -10.0],
		[0.0, 3.0, 0.0],
		[1.0, 3.0, 1.0],
		[1.0, -3.0, 1.0],
		[5.0, 0.0, 3.0],
		[3.0, 0.0, -5.0],
		[1.0, 3.0, 1.02],
		[1.02, 3.0, 1.0],
		[1.0, -3.0, -9.0]
    ]) {
        let a = test[0]
        let b = test[1]
        let c = test[2]
        let (rt1, rt2) = lFloat64.dlae2(a, b, c)

        let name = "a = ${a}, b = ${b}, c = ${c}"
		let a1 = a - rt1
		let c1 = c - rt1
		var det = a1 * c1 - b * b
        if (abs(det) > 1e-10) {
            throw Exception("First eigenvalue mismatch. ${name}. Det = ${det}")
        }

        let a2 = a - rt2
        let c2 = c - rt2
        det = a2 * c2 - b * b
        if (abs(det) > 1e-10) {
            throw Exception("Second eigenvalue mismatch. ${name}. Det = ${det}")
        }
    }
    return true
}