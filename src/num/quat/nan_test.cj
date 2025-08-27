package cjnum.num.quat

import std.unittest.*
import std.unittest.testmacro.*

private class NaNTest{
    let q:NqNumber
    let want:Bool
    public init(q:NqNumber,want:Bool) {
        this.q=q
        this.want=want
    }
}

private let nanTests:Array<NaNTest>=[
    NaNTest(nqNaN(), true),
    NaNTest(NqNumber(nan, nan, nan, nan), true),
    NaNTest(NqNumber(nan, 0.0, 0.0, 0.0), true),
    NaNTest(NqNumber(0.0, nan, 0.0, 0.0), true),
    NaNTest(NqNumber(0.0, 0.0, nan, 0.0), true),
    NaNTest(NqNumber(0.0, 0.0, 0.0, nan), true),
    NaNTest(NqNumber(inf, nan, nan, nan), false),
    NaNTest(NqNumber(nan, inf, nan, nan), false),
    NaNTest(NqNumber(nan, nan, inf, nan), false),
    NaNTest(NqNumber(nan, nan, nan, inf), false),
    NaNTest(NqNumber(-inf, nan, nan, nan), false),
    NaNTest(NqNumber(nan, -inf, nan, nan), false),
    NaNTest(NqNumber(nan, nan, -inf, nan), false),
    NaNTest(NqNumber(nan, nan, nan, -inf), false),
    NaNTest(NqNumber(0.0,0.0,0.0,0.0), false)
]

@Test
public func testIsNaN():Unit {
	for (test in nanTests) {
		let got = nqIsNaN(test.q)
		@Expect(got == test.want)
	}
}
