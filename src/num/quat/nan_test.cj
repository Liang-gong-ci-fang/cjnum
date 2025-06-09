package cjnum.num.quat

import std.math.*
import std.unittest.*
import std.unittest.testmacro.*

private class NaNTest{
    let q:Number
    let want:Bool
    public init(q:Number,want:Bool) {
        this.q=q
        this.want=want
    }
}

private let nanTests:Array<NaNTest>=[
    NaNTest(nqNaN(), true),
    NaNTest(Number(nan, nan, nan, nan), true),
    NaNTest(Number(nan, 0.0, 0.0, 0.0), true),
    NaNTest(Number(0.0, nan, 0.0, 0.0), true),
    NaNTest(Number(0.0, 0.0, nan, 0.0), true),
    NaNTest(Number(0.0, 0.0, 0.0, nan), true),
    NaNTest(Number(inf, nan, nan, nan), false),
    NaNTest(Number(nan, inf, nan, nan), false),
    NaNTest(Number(nan, nan, inf, nan), false),
    NaNTest(Number(nan, nan, nan, inf), false),
    NaNTest(Number(-inf, nan, nan, nan), false),
    NaNTest(Number(nan, -inf, nan, nan), false),
    NaNTest(Number(nan, nan, -inf, nan), false),
    NaNTest(Number(nan, nan, nan, -inf), false),
    NaNTest(Number(0.0,0.0,0.0,0.0), false)
]

@Test
public func testIsNaN():Unit {
	for (test in nanTests) {
		let got = nqIsNaN(test.q)
		@Expect(got == test.want)
	}
}
