package cjnum.floats.scalar

import std.math.*
import std.math.numeric.*
import std.unittest.*
import std.unittest.testmacro.*

const smallestNonzeroFloat64 = 4.9406564584124654e-324

class equalsRelativeTestStruct {
    public equalsRelativeTestStruct(
        let a: Float64,
        let b: Float64,
        var tol: Float64,
        let equal: Bool
    ) {}
}

@Test
public func testEqualsRelative() {
    let equalityTests = [
        equalsRelativeTestStruct(1000000.0, 1000001.0, 0.0, true),
        equalsRelativeTestStruct(1000001.0, 1000000.0, 0.0, true),
        equalsRelativeTestStruct(10000.0, 10001.0, 0.0, false),
        equalsRelativeTestStruct(10001.0, 10000.0, 0.0, false),
        equalsRelativeTestStruct(-1000000.0, -1000001.0, 0.0, true),
        equalsRelativeTestStruct(-1000001.0, -1000000.0, 0.0, true),
        equalsRelativeTestStruct(-10000.0, -10001.0, 0.0, false),
        equalsRelativeTestStruct(-10001.0, -10000.0, 0.0, false),
        equalsRelativeTestStruct(1.0000001, 1.0000002, 0.0, true),
        equalsRelativeTestStruct(1.0000002, 1.0000001, 0.0, true),
        equalsRelativeTestStruct(1.0002, 1.0001, 0.0, false),
        equalsRelativeTestStruct(1.0001, 1.0002, 0.0, false),
        equalsRelativeTestStruct(-1.000001, -1.000002, 0.0, true),
        equalsRelativeTestStruct(-1.000002, -1.000001, 0.0, true),
        equalsRelativeTestStruct(-1.0001, -1.0002, 0.0, false),
        equalsRelativeTestStruct(-1.0002, -1.0001, 0.0, false),
        equalsRelativeTestStruct(0.000000001000001, 0.000000001000002, 0.0, true),
        equalsRelativeTestStruct(0.000000001000002, 0.000000001000001, 0.0, true),
        equalsRelativeTestStruct(0.000000000001002, 0.000000000001001, 0.0, false),
        equalsRelativeTestStruct(0.000000000001001, 0.000000000001002, 0.0, false),
        equalsRelativeTestStruct(-0.000000001000001, -0.000000001000002, 0.0, true),
        equalsRelativeTestStruct(-0.000000001000002, -0.000000001000001, 0.0, true),
        equalsRelativeTestStruct(-0.000000000001002, -0.000000000001001, 0.0, false),
        equalsRelativeTestStruct(-0.000000000001001, -0.000000000001002, 0.0, false),
        equalsRelativeTestStruct(0.0, 0.0, 0.0, true),
        equalsRelativeTestStruct(0.0, -0.0, 0.0, true),
        equalsRelativeTestStruct(-0.0, -0.0, 0.0, true),
        equalsRelativeTestStruct(0.00000001, 0.0, 0.0, false),
        equalsRelativeTestStruct(0.0, 0.00000001, 0.0, false),
        equalsRelativeTestStruct(-0.00000001, 0.0, 0.0, false),
        equalsRelativeTestStruct(0.0, -0.00000001, 0.0, false),
        equalsRelativeTestStruct(0.0, 1e-310, 0.01, true),
        equalsRelativeTestStruct(1e-310, 0.0, 0.01, true),
        equalsRelativeTestStruct(1e-310, 0.0, 0.000001, false),
        equalsRelativeTestStruct(0.0, 1e-310, 0.000001, false),
        equalsRelativeTestStruct(0.0, -1e-310, 0.1, true),
        equalsRelativeTestStruct(-1e-310, 0.0, 0.1, true),
        equalsRelativeTestStruct(-1e-310, 0.0, 0.00000001, false),
        equalsRelativeTestStruct(0.0, -1e-310, 0.00000001, false),
        equalsRelativeTestStruct(Float64.Inf, Float64.Inf, 0.0, true),
        equalsRelativeTestStruct(-Float64.Inf, -Float64.Inf, 0.0, true),
        equalsRelativeTestStruct(-Float64.Inf, Float64.Inf, 0.0, false),
        equalsRelativeTestStruct(Float64.Inf, Float64.Max, 0.0, false),
        equalsRelativeTestStruct(-Float64.Inf, -Float64.Max, 0.0, false),
        equalsRelativeTestStruct(Float64.NaN, Float64.NaN, 0.0, false),
        equalsRelativeTestStruct(Float64.NaN, 0.0, 0.0, false),
        equalsRelativeTestStruct(-0.0, Float64.NaN, 0.0, false),
        equalsRelativeTestStruct(Float64.NaN, -0.0, 0.0, false),
        equalsRelativeTestStruct(0.0, Float64.NaN, 0.0, false),
        equalsRelativeTestStruct(Float64.NaN, Float64.Inf, 0.0, false),
        equalsRelativeTestStruct(Float64.Inf, Float64.NaN, 0.0, false),
        equalsRelativeTestStruct(Float64.NaN, -Float64.Inf, 0.0, false),
        equalsRelativeTestStruct(-Float64.Inf, Float64.NaN, 0.0, false),
        equalsRelativeTestStruct(Float64.NaN, Float64.Max, 0.0, false),
        equalsRelativeTestStruct(Float64.Max, Float64.NaN, 0.0, false),
        equalsRelativeTestStruct(Float64.NaN, -Float64.Max, 0.0, false),
        equalsRelativeTestStruct(-Float64.Max, Float64.NaN, 0.0, false),
        equalsRelativeTestStruct(Float64.NaN, smallestNonzeroFloat64, 0.0, false),
        equalsRelativeTestStruct(smallestNonzeroFloat64, Float64.NaN, 0.0, false),
        equalsRelativeTestStruct(Float64.NaN, -smallestNonzeroFloat64, 0.0, false),
        equalsRelativeTestStruct(-smallestNonzeroFloat64, Float64.NaN, 0.0, false),
        equalsRelativeTestStruct(1.000000001, -1.0, 0.0, false),
        equalsRelativeTestStruct(-1.0, 1.000000001, 0.0, false),
        equalsRelativeTestStruct(-1.000000001, 1.0, 0.0, false),
        equalsRelativeTestStruct(1.0, -1.000000001, 0.0, false),
        equalsRelativeTestStruct(10.0 * smallestNonzeroFloat64, 10.0 * (-smallestNonzeroFloat64), 0.0, true),
        equalsRelativeTestStruct(1e11 * smallestNonzeroFloat64, 1e11 * (-smallestNonzeroFloat64), 0.0, false),
        equalsRelativeTestStruct(smallestNonzeroFloat64, -smallestNonzeroFloat64, 0.0, true),
        equalsRelativeTestStruct(-smallestNonzeroFloat64, smallestNonzeroFloat64, 0.0, true),
        equalsRelativeTestStruct(smallestNonzeroFloat64, 0.0, 0.0, true),
        equalsRelativeTestStruct(0.0, smallestNonzeroFloat64, 0.0, true),
        equalsRelativeTestStruct(-smallestNonzeroFloat64, 0.0, 0.0, true),
        equalsRelativeTestStruct(0.0, -smallestNonzeroFloat64, 0.0, true),
        equalsRelativeTestStruct(0.000000001, -smallestNonzeroFloat64, 0.0, false),
        equalsRelativeTestStruct(0.000000001, smallestNonzeroFloat64, 0.0, false),
        equalsRelativeTestStruct(smallestNonzeroFloat64, 0.000000001, 0.0, false),
        equalsRelativeTestStruct(-smallestNonzeroFloat64, 0.000000001, 0.0, false)
    ]

    for (ts in equalityTests) {
        if (ts.tol == 0.0) {
            ts.tol = 1e-5
        }
        @Expect(equalWithinRel(ts.a, ts.b, ts.tol) == ts.equal)
    }
}

@Test
public func testEqualsULP() {
    @Expect(equalWithinULP(67329.242, nextAfterN(67329.242, Float64.Inf, 10), 10))
    @Expect(!equalWithinULP(67329.242, nextAfterN(67329.242, Float64.Inf, 5), 1))
    @Expect(!equalWithinULP(nextAfterN(67329.242, Float64.Inf, 5), 67329.242, 1))
    let f = nextAfterN(0.0, Float64.Inf, 2)
    @Expect(equalWithinULP(f, nextAfterN(f, -Float64.Inf, 5), 10))
    @Expect(equalWithinULP(67329.242, 67329.242, 10))
    @Expect(!equalWithinULP(1.0, Float64.NaN, 10))
}

func nextAfterN(x: Float64, y: Float64, n: Int64): Float64 {
    var xx = x
    for (i in 0..n) {
        xx = nextAfter(xx, y)
    }
    return xx
}

func nextAfter(x: Float64, y: Float64): Float64 {
    if (x.isNaN() || y.isNaN()) {
        return Float64.NaN
    }
    if (x == 0.0) {
        return copysign(Float64.fromBits(1), y)
    }
    if (x == Float64.Inf && y > x) {
        return x
    }
    if (x == -Float64.Inf && y < x) {
        return x
    }

    return if (x == y) {
        x
    } else if (x < y) {
        Float64.fromBits(x.toBits() + 1)
    } else {
        Float64.fromBits(x.toBits() - 1)
    }
}

class NaNWithTestStruct {
    NaNWithTestStruct(
        let payload: UInt64,
        var bits: UInt64
    ) {}
}

@Test
public func testNaNWith() {
    let nanWithTests = [
        NaNWithTestStruct(0, (0.0 / 0.0).toBits()),
        NaNWithTestStruct(0, Float64.NaN.toBits()),
        NaNWithTestStruct(1954, 0x7ff80000000007a2),
        NaNWithTestStruct(1, 9221120237041090561),
        NaNWithTestStruct(2, 9221120237041090562),
        NaNWithTestStruct(3, 9221120237041090563)
    ]

    for (test in nanWithTests) {
        let nan = naNWith(test.payload)
        @Expect(nan.isNaN())

        var bits = nan.toBits()

        const sign: UInt64 = 1 << 63
        bits = bits & (!sign)
        test.bits = test.bits & (!sign)

        @Expect(bits, test.bits)
    }
}

class NaNPayloadTestStruct {
    NaNPayloadTestStruct(
        let f: Float64,
        let payload: UInt64,
        let ok: Bool
    ) {}
}

@Test
public func testNaNPayload() {
    let nanPayloadTests = [
        NaNPayloadTestStruct(0.0 / 0.0, 0, true),
        NaNPayloadTestStruct(Float64.fromBits(Float64.NaN.toBits() | (1 << 63)), 0, true),
        NaNPayloadTestStruct(Float64.fromBits(Float64.NaN.toBits() & (!(1 << 63))), 0, true),
        NaNPayloadTestStruct(naNWith(1954), 1954, true),
        NaNPayloadTestStruct(copysign(0.0, -1.0), 0, false),
        NaNPayloadTestStruct(0.0, 0, false),
        NaNPayloadTestStruct(-Float64.Inf, 0, false),
        NaNPayloadTestStruct(Float64.Inf, 0, false),
        NaNPayloadTestStruct(Float64.fromBits(0x7ff0000000000001), 0, false)
    ]

    for (test in nanPayloadTests) {
        let (payload, ok) = naNPayload(test.f)
        @Expect(payload, test.payload)
        @Expect(ok, test.ok)
    }
}

func copysign(x: Float64, y: Float64): Float64 {
    let sign = if (y > 0.0) {
        1.0
    } else {
        -1.0
    }
    return sign * x
}

class roundTestStruct {
    roundTestStruct(
        let x: Float64,
        let prec: Int64,
        let want: Float64
    ) {}
}

@Test
public func testRound() {
    let roundTests = [
        roundTestStruct(0.0, 1, 0.0),
        roundTestStruct(Float64.Inf, 1, Float64.Inf),
        roundTestStruct(-Float64.Inf, 1, -Float64.Inf),
        roundTestStruct(Float64.NaN, 1, Float64.NaN),
        roundTestStruct(-0.0, 1, 0.0),
        roundTestStruct(Float64.Max / 2.0, 1, Float64.Max / 2.0),
        roundTestStruct(Decimal(2).powWithPrecision(64, 0).toFloat64(), 1, Decimal(2).powWithPrecision(64, 0).toFloat64(
        )),
        roundTestStruct(454.4445, 3, 454.445),
        roundTestStruct(454.44445, 4, 454.4445),
        roundTestStruct(0.42499, 4, 0.425),
        roundTestStruct(0.42599, 4, 0.426),
        roundTestStruct(0.424999999999993, 2, 0.42),
        roundTestStruct(0.425, 2, 0.43),
        roundTestStruct(0.425000000000001, 2, 0.43),
        roundTestStruct(123.4244999999999, 3, 123.424),
        roundTestStruct(123.4245, 3, 123.425),
        roundTestStruct(123.4245000000001, 3, 123.425),
        roundTestStruct(-0.49999999999999994, 0, 0.0),
        roundTestStruct(0.49999999999999994, 0, 0.0),
        roundTestStruct(454.45, 0, 454.0),
        roundTestStruct(454.45, 1, 454.5),
        roundTestStruct(454.45, 2, 454.45),
        roundTestStruct(454.45, 3, 454.45),
        roundTestStruct(454.445, 0, 454.0),
        roundTestStruct(454.445, 1, 454.4),
        roundTestStruct(454.445, 2, 454.45),
        roundTestStruct(454.445, 3, 454.445),
        roundTestStruct(454.445, 4, 454.445),
        roundTestStruct(454.55, 0, 455.0),
        roundTestStruct(454.55, 1, 454.6),
        roundTestStruct(454.55, 2, 454.55),
        roundTestStruct(454.55, 3, 454.55),
        roundTestStruct(454.455, 0, 454.0),
        roundTestStruct(454.455, 1, 454.5),
        roundTestStruct(454.455, 2, 454.46),
        roundTestStruct(454.455, 3, 454.455),
        roundTestStruct(454.455, 4, 454.455),
        roundTestStruct(Float64.Inf, -1, Float64.Inf),
        roundTestStruct(-Float64.Inf, -1, -Float64.Inf),
        roundTestStruct(Float64.NaN, -1, Float64.NaN),
        roundTestStruct(-0.0, -1, 0.0),
        roundTestStruct(454.45, -1, 450.0),
        roundTestStruct(454.45, -2, 500.0),
        roundTestStruct(500.0, -3, 1000.0),
        roundTestStruct(500.0, -4, 0.0),
        roundTestStruct(1500.0, -3, 2000.0),
        roundTestStruct(1500.0, -4, 0.0)
    ]

    for (test in roundTests) {
        for (sign in [1.0, -1.0]) {
            let got = round(sign * test.x, test.prec)
            var want = sign * test.want
            if (want == 0.0) {
                want = 0.0
            }
            @Expect(!((got != want || got >= 0.0 != want >= 0.0) && !(got.isNaN() && want.isNaN())))
        }
    }
}

@Test
public func testRoundEven() {
    let roundEvenTests = [
        roundTestStruct(0.0, 1, 0.0),
        roundTestStruct(Float64.Inf, 1, Float64.Inf),
        roundTestStruct(-Float64.Inf, 1, -Float64.Inf),
        roundTestStruct(Float64.NaN, 1, Float64.NaN),
        roundTestStruct(-0.0, 1, 0.0),
        roundTestStruct(Float64.Max / 2.0, 1, Float64.Max / 2.0),
        roundTestStruct(Decimal(2).powWithPrecision(64, 0).toFloat64(), 1, Decimal(2).powWithPrecision(64, 0).toFloat64(
        )),
        roundTestStruct(454.4445, 3, 454.444),
        roundTestStruct(454.44445, 4, 454.4444),
        roundTestStruct(0.42499, 4, 0.425),
        roundTestStruct(0.42599, 4, 0.426),
        roundTestStruct(0.424999999999993, 2, 0.42),
        roundTestStruct(0.425, 2, 0.42),
        roundTestStruct(0.425000000000001, 2, 0.43),
        roundTestStruct(123.4244999999999, 3, 123.424),
        roundTestStruct(123.4245, 3, 123.424),
        roundTestStruct(123.4245000000001, 3, 123.425),
        roundTestStruct(-0.49999999999999994, 0, 0.0),
        roundTestStruct(0.49999999999999994, 0, 0.0),
        roundTestStruct(454.45, 0, 454.0),
        roundTestStruct(454.45, 1, 454.4),
        roundTestStruct(454.45, 2, 454.45),
        roundTestStruct(454.45, 3, 454.45),
        roundTestStruct(454.445, 0, 454.0),
        roundTestStruct(454.445, 1, 454.4),
        roundTestStruct(454.445, 2, 454.44),
        roundTestStruct(454.445, 3, 454.445),
        roundTestStruct(454.445, 4, 454.445),
        roundTestStruct(454.55, 0, 455.0),
        roundTestStruct(454.55, 1, 454.6),
        roundTestStruct(454.55, 2, 454.55),
        roundTestStruct(454.55, 3, 454.55),
        roundTestStruct(454.455, 0, 454.0),
        roundTestStruct(454.455, 1, 454.5),
        roundTestStruct(454.455, 2, 454.46),
        roundTestStruct(454.455, 3, 454.455),
        roundTestStruct(454.455, 4, 454.455),
        roundTestStruct(Float64.Inf, -1, Float64.Inf),
        roundTestStruct(-Float64.Inf, -1, -Float64.Inf),
        roundTestStruct(Float64.NaN, -1, Float64.NaN),
        roundTestStruct(-0.0, -1, 0.0),
        roundTestStruct(454.45, -1, 450.0),
        roundTestStruct(454.45, -2, 500.0),
        roundTestStruct(500.0, -3, 0.0),
        roundTestStruct(500.0, -4, 0.0),
        roundTestStruct(1500.0, -3, 2000.0),
        roundTestStruct(1500.0, -4, 0.0)
    ]

    for (test in roundEvenTests) {
        for (sign in [1.0, -1.0]) {
            let got = roundEven(sign * test.x, test.prec)
            var want = sign * test.want
            if (want == 0.0) {
                want = 0.0
            }
            @Expect(!((got != want || got >= 0.0 != want >= 0.0) && !(got.isNaN() && want.isNaN())))
        }
    }
}

class sameTestStruct {
    sameTestStruct(
        let a: Float64,
        let b: Float64,
        let want: Bool
    ) {}
}

@Test
public func sameTest() {
    let sameTests = [
        sameTestStruct(0.0, 0.0, true),
        sameTestStruct(1.0, 1.0, true),
        sameTestStruct(-1.0, 1.0, false),
        sameTestStruct(0.0, 1.0, false),
        sameTestStruct(1.0, 0.0, false),
        sameTestStruct(-1.0, 1.0, false),
        sameTestStruct(Float64.NaN, Float64.NaN, true),
        sameTestStruct(1.0, Float64.NaN, false),
        sameTestStruct(Float64.Inf, Float64.NaN, false),
        sameTestStruct(Float64.NaN, Float64.Inf, false),
        sameTestStruct(Float64.NaN, 1.0, false),
        sameTestStruct(Float64.Inf, Float64.Inf, true),
        sameTestStruct(-Float64.Inf, Float64.NaN, false),
        sameTestStruct(Float64.Inf, -Float64.Inf, false)
    ]

    for (test in sameTests) {
        let got = same(test.a, test.b)
        @Expect(got, test.want)
    }
}
