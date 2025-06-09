package cjnum.num.quat

import std.math.*
import std.unittest.*
import std.unittest.testmacro.*

private class ExpTest {
    let q: Number
    let want: Number
    public init(q: Number, want: Number) {
        this.q = q
        this.want = want
    }
}

private let expTests: Array<ExpTest> = [
    ExpTest(
        Number(0.0, 0.0, 0.0, 0.0),
        Number(1.0, 0.0, 0.0, 0.0),
    ),
    ExpTest(
        Number(1.0, 1.0, 1.0, 1.0),
        Number(-0.43643792124786496, 1.549040352371697, 1.549040352371697, 1.549040352371697),
    ),
    ExpTest(
        Number(1.0, 0.0, 1.0, 1.0),
        Number(0.42389891174348104, 0.0, 1.8986002490721081, 1.8986002490721081),
    ),
    ExpTest(
        Number(1.0, 0.0, 0.0, 1.0),
        Number(1.4686939399158851, 0.0, 0.0, 2.2873552871788423),
    ),
    ExpTest(
        Number(0.0, 1.0, 1.0, 1.0),
        Number(-0.16055653857469052, 0.569860099182514, 0.569860099182514, 0.569860099182514),
    ),
    ExpTest(
        Number(0.0, 0.0, 1.0, 1.0),
        Number(0.15594369476537437, 0.0, 0.6984559986366083, 0.6984559986366083),
    ),
    ExpTest(
        Number(0.0, 0.0, 0.0, 1.0),
        Number(0.5403023058681398, 0.0, 0.0, 0.8414709848078965),
    )
]

@Test
public func testExp(): Unit {
    let tol = 1e-14
    for (test in expTests) {
        let got = nqExp(test.q)
        @Expect(nqEqualApprox(got, test.want, tol))
    }
}

private class LogTest {
    let q: Number
    let want: Number
    public init(q: Number, want: Number) {
        this.q = q
        this.want = want
    }
}

private let logTests: Array<LogTest> = [
    LogTest(
        Number(0.0, 0.0, 0.0, 0.0),
        Number(-inf, 0.0, 0.0, 0.0),
    ),
    LogTest(
        Number(1.0, 1.0, 1.0, 1.0),
        Number(0.6931471805599453, 0.6045997880780728, 0.6045997880780728, 0.6045997880780728),
    ),
    LogTest(
        Number(1.0, 0.0, 1.0, 1.0),
        Number(0.5493061443340548, 0.0, 0.6755108588560398, 0.6755108588560398),
    ),
    LogTest(
        Number(1.0, 0.0, 0.0, 1.0),
        Number(0.3465735902799727, 0.0, 0.0, 0.7853981633974484),
    ),
    LogTest(
        Number(0.0, 1.0, 1.0, 1.0),
        Number(0.5493061443340548, 0.906899682117109, 0.906899682117109, 0.906899682117109),
    ),
    LogTest(
        Number(0.0, 0.0, 1.0, 1.0),
        Number(0.3465735902799727, 0.0, 1.1107207345395915, 1.1107207345395915),
    ),
    LogTest(
        Number(0.0, 0.0, 0.0, 1.0),
        Number(0.0, 0.0, 0.0, 1.5707963267948966),
    )
]

@Test
public func testLog(): Unit {
    let tol = 1e-14
    for (test in logTests) {
        let got = nqLog(test.q)
        @Expect(nqEqualApprox(got, test.want, tol))
    }
}

public class PowTest {
    let q: Number
    let r: Number
    let want: Number
    public init(q: Number, r: Number, want: Number) {
        this.q = q
        this.r = r
        this.want = want
    }
}

private let powTests: Array<PowTest> = [
    PowTest(Number(0.0, 0.0, 0.0, 0.0), Number(0.0, 0.0, 0.0, 0.0), Number(1.0, 0.0, 0.0, 0.0)),
    // Expected velues below are from pyquaternion.
    // pyquaternion does not support quaternion powers.
    // TODO(kortschak): Add non-real r cases.
    PowTest(Number(1.0, 1.0, 1.0, 1.0), Number(2.0, 0.0, 0.0, 0.0), Number(-2.0, 2.0, 2.0, 2.0)),
    PowTest(Number(1.0, 0.0, 1.0, 1.0), Number(2.0, 0.0, 0.0, 0.0), Number(-1.0, 0.0, 2.0, 2.0)),
    PowTest(Number(1.0, 0.0, 0.0, 1.0), Number(2.0, 0.0, 0.0, 0.0), Number(0.0, 0.0, 0.0, 2.0)),
    PowTest(Number(0.0, 1.0, 1.0, 1.0), Number(2.0, 0.0, 0.0, 0.0), Number(-3.0, 0.0, 0.0, 0.0)),
    PowTest(Number(0.0, 0.0, 1.0, 1.0), Number(2.0, 0.0, 0.0, 0.0), Number(-2.0, 0.0, 0.0, 0.0)),
    PowTest(Number(0.0, 0.0, 0.0, 1.0), Number(2.0, 0.0, 0.0, 0.0), Number(-1.0, 0.0, 0.0, 0.0)),
    PowTest(Number(1.0, 1.0, 1.0, 1.0), Number(Float64.PI, 0.0, 0.0, 0.0),
        Number(-8.728144138959564, -0.7527136547040768, -0.7527136547040768, -0.7527136547040768)),
    PowTest(Number(1.0, 0.0, 1.0, 1.0), Number(Float64.PI, 0.0, 0.0, 0.0),
        Number(-5.561182514695044, 0.0, 0.5556661490713818, 0.5556661490713818)),
    PowTest(Number(1.0, 0.0, 0.0, 1.0), Number(Float64.PI, 0.0, 0.0, 0.0),
        Number(-2.320735561810013, 0.0, 0.0, 1.8544983901925216)),
    PowTest(Number(0.0, 1.0, 1.0, 1.0), Number(Float64.PI, 0.0, 0.0, 0.0),
        Number(1.2388947209955585, -3.162774128856231, -3.162774128856231, -3.162774128856231)),
    PowTest(Number(0.0, 0.0, 1.0, 1.0), Number(Float64.PI, 0.0, 0.0, 0.0),
        Number(0.6552860151073727, 0.0, -2.0488506614051922, -2.0488506614051922)),
    PowTest(Number(0.0, 0.0, 0.0, 1.0), Number(Float64.PI, 0.0, 0.0, 0.0),
        Number(0.22058404074969779, 0.0, 0.0, -0.9753679720836315)),
    PowTest(Number(1.0, 1.0, 1.0, 1.0), Number(3.0, 0.0, 0.0, 0.0), Number(-8.0, 0.0, 0.0, 0.0)),
    PowTest(Number(1.0, 0.0, 1.0, 1.0), Number(3.0, 0.0, 0.0, 0.0), Number(-5.0, 0.0, 1.0, 1.0)),
    PowTest(Number(1.0, 0.0, 0.0, 1.0), Number(3.0, 0.0, 0.0, 0.0), Number(-2.0, 0.0, 0.0, 2.0)),
    PowTest(Number(0.0, 1.0, 1.0, 1.0), Number(3.0, 0.0, 0.0, 0.0), Number(0.0, -3.0, -3.0, -3.0)),
    PowTest(Number(0.0, 0.0, 1.0, 1.0), Number(3.0, 0.0, 0.0, 0.0), Number(0.0, 0.0, -2.0, -2.0)),
    PowTest(Number(0.0, 0.0, 0.0, 1.0), Number(3.0, 0.0, 0.0, 0.0), Number(0.0, 0.0, 0.0, -1.0),)
]

@Test
public func testPow(): Unit {
    let tol = 1e-14
    for (test in powTests) {
        let got = nqPow(test.q, test.r)
        @Expect(nqEqualApprox(got, test.want, tol))
    }
}

private class PowRealTest {
    let q: Number
    let r: Float64
    let want: Number
    public init(q: Number, r: Float64, want: Number) {
        this.q = q
        this.r = r
        this.want = want
    }
}

private let powRealTests: Array<PowRealTest> = [
    PowRealTest(Number(0.0, 0.0, 0.0, 0.0), 0.0, Number(1.0, 0.0, 0.0, 0.0)),
    PowRealTest(Number(1.0, 1.0, 1.0, 1.0), 2.0, Number(-2.0, 2.0, 2.0, 2.0)),
    PowRealTest(Number(1.0, 0.0, 1.0, 1.0), 2.0, Number(-1.0, 0.0, 2.0, 2.0)),
    PowRealTest(Number(1.0, 0.0, 0.0, 1.0), 2.0, Number(0.0, 0.0, 0.0, 2.0)),
    PowRealTest(Number(0.0, 1.0, 1.0, 1.0), 2.0, Number(-3.0, 0.0, 0.0, 0.0)),
    PowRealTest(Number(0.0, 0.0, 1.0, 1.0), 2.0, Number(-2.0, 0.0, 0.0, 0.0)),
    PowRealTest(Number(0.0, 0.0, 0.0, 1.0), 2.0, Number(-1.0, 0.0, 0.0, 0.0)),
    PowRealTest(Number(1.0, 1.0, 1.0, 1.0), Float64.PI,
        Number(-8.728144138959564, -0.7527136547040768, -0.7527136547040768, -0.7527136547040768)),
    PowRealTest(Number(1.0, 0.0, 1.0, 1.0), Float64.PI,
        Number(-5.561182514695044, 0.0, 0.5556661490713818, 0.5556661490713818)),
    PowRealTest(Number(1.0, 0.0, 0.0, 1.0), Float64.PI, Number(-2.320735561810013, 0.0, 0.0, 1.8544983901925216)),
    PowRealTest(Number(0.0, 1.0, 1.0, 1.0), Float64.PI,
        Number(1.2388947209955585, -3.162774128856231, -3.162774128856231, -3.162774128856231)),
    PowRealTest(Number(0.0, 0.0, 1.0, 1.0), Float64.PI,
        Number(0.6552860151073727, 0.0, -2.0488506614051922, -2.0488506614051922)),
    PowRealTest(Number(0.0, 0.0, 0.0, 1.0), Float64.PI, Number(0.22058404074969779, 0.0, 0.0, -0.9753679720836315)),
    PowRealTest(Number(1.0, 1.0, 1.0, 1.0), 3.0, Number(-8.0, 0.0, 0.0, 0.0)),
    PowRealTest(Number(1.0, 0.0, 1.0, 1.0), 3.0, Number(-5.0, 0.0, 1.0, 1.0)),
    PowRealTest(Number(1.0, 0.0, 0.0, 1.0), 3.0, Number(-2.0, 0.0, 0.0, 2.0)),
    PowRealTest(Number(0.0, 1.0, 1.0, 1.0), 3.0, Number(0.0, -3.0, -3.0, -3.0)),
    PowRealTest(Number(0.0, 0.0, 1.0, 1.0), 3.0, Number(0.0, 0.0, -2.0, -2.0)),
    PowRealTest(Number(0.0, 0.0, 0.0, 1.0), 3.0, Number(0.0, 0.0, 0.0, -1.0))
]

@Test
public func testRealPow(): Unit {
    let tol = 1e-14
    for (test in powRealTests) {
        let got = nqPowReal(test.q, test.r)
        @Expect(nqEqualApprox(got, test.want, tol))
    }
}

private class SqrtTest {
    let q: Number
    let want: Number
    public init(q: Number, want: Number) {
        this.q = q
        this.want = want
    }
}

private let sqrtTests: Array<SqrtTest> = [
    SqrtTest(Number(0.0, 0.0, 0.0, 0.0), Number(0.0, 0.0, 0.0, 0.0)),
    SqrtTest(Number(1.0, 1.0, 1.0, 1.0),
        Number(1.2247448713915892, 0.4082482904638631, 0.4082482904638631, 0.4082482904638631)),
    SqrtTest(Number(1.0, 0.0, 1.0, 1.0), Number(1.1687708944803676, 0.0, 0.42779983858367593, 0.42779983858367593)),
    SqrtTest(Number(1.0, 0.0, 0.0, 1.0), Number(1.0986841134678098, 0.0, 0.0, 0.45508986056222733)),
    SqrtTest(Number(0.0, 1.0, 1.0, 1.0),
        Number(0.9306048591020996, 0.5372849659117709, 0.5372849659117709, 0.5372849659117709)),
    SqrtTest(Number(0.0, 0.0, 1.0, 1.0), Number(0.8408964152537146, 0.0, 0.5946035575013604, 0.5946035575013604)),
    SqrtTest(Number(0.0, 0.0, 0.0, 1.0), Number(0.7071067811865476, 0.0, 0.0, 0.7071067811865475))
]

@Test
public func testSqrt(): Unit {
    let tol = 1e-14
    for (test in sqrtTests) {
        let got = nqSqrt(test.q)
        @Expect(nqEqualApprox(got, test.want, tol))
    }
}
