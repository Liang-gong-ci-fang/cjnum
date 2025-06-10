package cjnum.lapack.testlapack

import cjnum.lapack.*
import std.math.*
import std.random.*

public func dlaev2Test(lFloat64: LFloat64): Bool {
    let rnd = Random(1)
    for (_ in 0..100) {
        let a = rnd.nextGaussianFloat64()
        let b = rnd.nextGaussianFloat64()
        let c = rnd.nextGaussianFloat64()

        let (rt1, rt2, cs1, sn1) = lFloat64.dlaev2(a, b, c)
        let tmp = mul2by2([[cs1, sn1], [-sn1, cs1]], [[a, b], [b, c]])
        let ans = mul2by2(tmp, [[cs1, -sn1], [sn1, cs1]])
        if (abs(ans[0][0] - rt1) > 1e-14) {
            throw Exception("Largest eigenvalue mismatch. Returned ${rt1}, mul ${ans[0][0]}")
        }
        if (abs(ans[1][0]) > 1e-14 || abs(ans[0][1]) > 1e-14) {
            throw Exception("Non-zero off diagonal. ans[1][0] = ${ans[1][0]}, ans[0][1] = ${ans[0][1]}")
        }
        if (abs(ans[1][1] - rt2) > 1e-14) {
            throw Exception("Smallest eigenvalue mismatch. Returned ${rt2}, mul ${ans[1][1]}")
        }
    }
    return true
}

func mul2by2(a: Array<Array<Float64>>, b: Array<Array<Float64>>): Array<Array<Float64>> {
    var c: Array<Array<Float64>> = [[0.0, 0.0], [0.0, 0.0]]
    c[0][0] = a[0][0] * b[0][0] + a[0][1] * b[1][0]
    c[0][1] = a[0][0] * b[0][1] + a[0][1] * b[1][1]
    c[1][0] = a[1][0] * b[0][0] + a[1][1] * b[1][0]
    c[1][1] = a[1][0] * b[0][1] + a[1][1] * b[1][1]
    return c
}
