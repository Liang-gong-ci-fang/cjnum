package cjnum.num.quat

import std.math.*
import std.unittest.*
import std.unittest.testmacro.*

private class AbsTest {
    let q: Number
    let want: Float64
    public init(q: Number, want: Float64) {
        this.q = q
        this.want = want
    }
}

private let absTests: Array<AbsTest> = [
    AbsTest(Number(0.0, 0.0, 0.0, 0.0), 0.0),
    AbsTest(nqNaN(), nan),
    AbsTest(nqInf(), inf),
    AbsTest(Number(1.0, 1.0, 1.0, 1.0), 2.0),
    AbsTest(Number(-1.0, 1.0, -1.0, 1.0), 2.0),
    AbsTest(Number(1.0, 2.0, 3.0, 4.0), sqrt(1.0 + 4.0 + 9.0 + 16.0)),
    AbsTest(Number(-1.0, -2.0, -3.0, -4.0), sqrt(1.0 + 4.0 + 9.0 + 16.0))
]

@Test
public func testAbs(): Unit {
    for (test in absTests) {
        let got = nqAbs(test.q)
        if (test.want.isNaN()) {
            @Expect(got.isNaN())
            continue
        }
        @Expect(got==test.want)
    }
}
