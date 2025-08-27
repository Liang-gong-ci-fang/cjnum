package cjnum.cmplxs.cscalar

import cjnum.complex.*
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
        equalsRelativeTestStruct(0.0, 1e-310, 1.0, true),
        equalsRelativeTestStruct(1e-310, 0.0, 1.0, true),
        equalsRelativeTestStruct(1e-310, 0.0, 0.000001, false),
        equalsRelativeTestStruct(0.0, 1e-310, 0.000001, false),
        equalsRelativeTestStruct(0.0, -1e-310, 1.0, true),
        equalsRelativeTestStruct(-1e-310, 0.0, 1.0, true),
        equalsRelativeTestStruct(-1e-310, 0.0, 0.00000001, false),
        equalsRelativeTestStruct(0.0, -1e-310, 0.00000001, false),
        equalsRelativeTestStruct(Float64.Inf, Float64.Inf, 0.0, true),
        equalsRelativeTestStruct(Float64.Inf, Float64.Max, 0.0, false),
        equalsRelativeTestStruct(-Float64.Inf, -Float64.Max, 0.0, false),
        equalsRelativeTestStruct(Float64.NaN, Float64.NaN, 0.0, false),
        equalsRelativeTestStruct(Float64.NaN, 0.0, 0.0, false),
        equalsRelativeTestStruct(-0.0, Float64.NaN, 0.0, false),
        equalsRelativeTestStruct(Float64.NaN, -0.0, 0.0, false),
        equalsRelativeTestStruct(0.0, Float64.NaN, 0.0, false),
        equalsRelativeTestStruct(Float64.NaN, Float64.Inf, 0.0, false),
        equalsRelativeTestStruct(Float64.Inf, Float64.NaN, 0.0, false),
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
        let complexToCompare = [
            (Complex128(ts.a, 0.0), Complex128(ts.b, 0.0)),
            (Complex128(0.0, ts.a), Complex128(0.0, ts.b)),
            (Complex128(ts.a, ts.a), Complex128(ts.b, ts.b))
        ]
        for (it in complexToCompare){
        @Expect(equalWithinRel(it[0], it[1], ts.tol) == ts.equal)
        }
    }
}

class roundTestStruct {
    roundTestStruct(
        let x: Complex128,
        let prec: Int64,
        let want: Complex128
    ) {}
}

@Test
public func testRound() {
    let roundTests = [
        roundTestStruct(Complex128(0.0), 1, Complex128(0.0)),
        roundTestStruct(Complex128(Float64.Inf, Float64.Inf), 1, Complex128(Float64.Inf, Float64.Inf)),
        roundTestStruct(Complex128(Float64.NaN, Float64.NaN), 1, Complex128(Float64.NaN, Float64.NaN)),
        roundTestStruct(Complex128(0.0), 1, Complex128(0.0)),
        roundTestStruct(Complex128(Float64.Max/2.0), 1, Complex128(Float64.Max/2.0)),
        roundTestStruct(Complex128(1.844674e19), 1, Complex128(1.844674e19)),
        roundTestStruct(Complex128(454.4445), 3, Complex128(454.444)),
        roundTestStruct(Complex128(454.44445), 4, Complex128(454.4444)),
        roundTestStruct(Complex128(0.42499), 4, Complex128(0.425)),
        roundTestStruct(Complex128(0.42599), 4, Complex128(0.426)),
        roundTestStruct(Complex128(0.424999999999993), 2, Complex128(0.42)),
        roundTestStruct(Complex128(0.425), 2, Complex128(0.42)),
        roundTestStruct(Complex128(0.425000000000001), 2, Complex128(0.43)),
        roundTestStruct(Complex128(123.4244999999999), 3, Complex128(123.424)),
        roundTestStruct(Complex128(123.4245), 3, Complex128(123.424)),
        roundTestStruct(Complex128(123.4245000000001), 3, Complex128(123.425)),
        roundTestStruct(Complex128(454.45), 0, Complex128(454.0)),
        roundTestStruct(Complex128(454.45), 1, Complex128(454.4)),
        roundTestStruct(Complex128(454.45), 2, Complex128(454.45)),
        roundTestStruct(Complex128(454.45), 3, Complex128(454.45)),
        roundTestStruct(Complex128(454.445), 0, Complex128(454.0)),
        roundTestStruct(Complex128(454.445), 1, Complex128(454.4)),
        roundTestStruct(Complex128(454.445), 2, Complex128(454.44)),
        roundTestStruct(Complex128(454.445), 3, Complex128(454.445)),
        roundTestStruct(Complex128(454.445), 4, Complex128(454.445)),
        roundTestStruct(Complex128(454.55), 0, Complex128(455.0)),
        roundTestStruct(Complex128(454.55), 1, Complex128(454.6)),
        roundTestStruct(Complex128(454.55), 2, Complex128(454.55)),
        roundTestStruct(Complex128(454.55), 3, Complex128(454.55)),
        roundTestStruct(Complex128(454.455), 0, Complex128(454.0)),
        roundTestStruct(Complex128(454.455), 1, Complex128(454.5)),
        roundTestStruct(Complex128(454.455), 2, Complex128(454.46)),
        roundTestStruct(Complex128(454.455), 3, Complex128(454.455)),
        roundTestStruct(Complex128(454.455), 4, Complex128(454.455)),
        roundTestStruct(Complex128(454.455), -1, Complex128(450.0)),
        roundTestStruct(Complex128(454.455), -2, Complex128(500.0)),
        roundTestStruct(Complex128(500.0), -3, Complex128(0.0)),
        roundTestStruct(Complex128(500.0), -4, Complex128(0.0)),
        roundTestStruct(Complex128(1500.0), -3, Complex128(2000.0)),
        roundTestStruct(Complex128(1500.0), -4, Complex128(0.0))
    ]

    for(sign in roundTests) {
        let result = roundEven(sign.x, sign.prec)
        @Expect(same(result, sign.want))
    }
}

class sameTestStruct {
    sameTestStruct(
        let a: Complex128,
        let b: Complex128,
        let want: Bool
    ){}
}

@Test
public func testSame() {
    let sameTests = [
        sameTestStruct(Complex128(0.0), Complex128(0.0), true),
        sameTestStruct(Complex128(1.0), Complex128(1.0), true),
        sameTestStruct(Complex128(-1.0), Complex128(1.0), false),
        sameTestStruct(Complex128(0.0), Complex128(1.0), false),
        sameTestStruct(Complex128(1.0), Complex128(0.0), false),
        sameTestStruct(Complex128(-1.0), Complex128(1.0), false),
        sameTestStruct(Complex128(Float64.NaN, Float64.NaN), Complex128(Float64.NaN, Float64.NaN), true),
        sameTestStruct(Complex128(1.0), Complex128(Float64.NaN, Float64.NaN), false),
        sameTestStruct(Complex128(Float64.Inf, Float64.Inf), Complex128(Float64.NaN, Float64.NaN), false),
        sameTestStruct(Complex128(Float64.NaN, Float64.NaN), Complex128(0.0), false),
        sameTestStruct(Complex128(Float64.NaN, Float64.NaN), Complex128(1.0), false),
        sameTestStruct(Complex128(Float64.Inf, Float64.Inf), Complex128(Float64.Inf, Float64.Inf), true)
    ]
    
    for (test in sameTests) {
        let result = same(test.a, test.b)
        @Expect(result == test.want)
    }
}