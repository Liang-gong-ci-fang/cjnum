package cjnum.lapack.testlapack

import std.math.*
import cjnum.floats.*
import cjnum.lapack.*

class drsclTestStruct {
    drsclTestStruct(
        let x!: Array<Float64>,
        let a!: Float64
    ) { }
}

public func drsclTest(lFloat64: LFloat64): Bool {
    let tests = [
        drsclTestStruct(
            x: [1.0, 2.0, 3.0, 4.0, 5.0],
            a: 4.0,
        ),
        drsclTestStruct(
            x: [1.0, 2.0, 3.0, 4.0, 5.0],
            a: Float64.Max
        ),
        drsclTestStruct(
            x: [1.0, 2.0, 3.0, 4.0, 5.0],
            a: 1e-307
        )
    ]
    for (i in 0..tests.size) {
        let test = tests[i]
        let xcopy = test.x.clone()

        lFloat64.drscl(test.x.size, test.a, xcopy, 1)
        if (fEqual(xcopy, test.x)) {
            throw Exception("x unchanged during call to drscl. a = ${test.a}, x = ${test.x}.")
        }
        lFloat64.drscl(test.x.size, 1.0 / test.a, xcopy, 1)
        if (!fEqualApprox(xcopy, test.x, 1e-14)) {
            throw Exception("x not equal after scaling and unscaling. a = ${test.a}, x = ${test.x}.")
        }
    }
    return true
}