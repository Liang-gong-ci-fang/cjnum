package cjnum.num.quat

import std.math.*
import std.unittest.*
import std.unittest.testmacro.*

private class InfTest{
    let q:NqNumber
    let want:Bool
    public init(q:NqNumber,want:Bool){
        this.q=q
        this.want=want
    }
}

private let infTests:Array<InfTest>=[
    InfTest(nqInf(),true),
    InfTest(NqNumber(inf, inf, inf, inf), true),
    InfTest(NqNumber(-inf, -inf, -inf, -inf), true),
    InfTest(NqNumber(inf, nan, nan, nan), true),
    InfTest(NqNumber(nan, inf, nan, nan), true),
    InfTest(NqNumber(nan, nan, inf, nan), true),
    InfTest(NqNumber(nan, nan, nan, inf), true),
    InfTest(NqNumber(-inf, nan, nan, nan), true),
    InfTest(NqNumber(nan, -inf, nan, nan), true),
    InfTest(NqNumber(nan, nan, -inf, nan), true),
    InfTest(NqNumber(nan, nan, nan, -inf), true),
    InfTest(NqNumber(inf, 0.0, 0.0, 0.0), true),
    InfTest(NqNumber(0.0, inf, 0.0, 0.0), true),
    InfTest(NqNumber(0.0, 0.0, inf, 0.0), true),
    InfTest(NqNumber(0.0, 0.0, 0.0, inf), true),
    InfTest(NqNumber(-inf, 0.0, 0.0, 0.0), true),
    InfTest(NqNumber(0.0, -inf, 0.0, 0.0), true),
    InfTest(NqNumber(0.0, 0.0, -inf, 0.0), true),
    InfTest(NqNumber(0.0, 0.0, 0.0, -inf), true),
    InfTest(NqNumber(0.0,0.0,0.0,0.0), false)
]

@Test
public func testIsInf():Unit {
	for (test in infTests) {
		let got = nqIsInf(test.q)
        @Expect(got == test.want)
	}
}