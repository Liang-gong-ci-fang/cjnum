package cjnum.num.dualquat

import cjnum.util.*
import cjnum.num.quat.*
import cjnum.floats.scalar.*
import std.math.*
import std.unittest.*
import std.unittest.testmacro.*

// TODO var formatTests = []struct {

// TODO func TestFormat(t *testing.T) {

// First derivatives:

private func dExp(x: NqNumber): NqNumber {
    return nqExp(x)
}

private func dLog(x: NqNumber): NqNumber {
    if (x.nqEqual(zeroQuat)) {
        return nqInf()
    } else if (nqIsInf(x)) {
        return zeroQuat
    }
    return nqInv(x)
}

private func dSqrt(x: NqNumber): NqNumber {
    return nqScale(0.5, nqInv(nqSqrt(x)))
}

private func dInv(x: NqNumber): NqNumber {
    return nqScale(-1.0, nqInv(nqMul(x, x)))
}

private let negZero = copysign(0.0, -1.0)
private let oneReal = NqNumber(1.0, 0.0, 0.0, 0.0)
private let negZeroQuat = nqScale(-1.0, NqNumber(0.0, 0.0, 0.0, 0.0))
private let one = NqNumber(1.0, 1.0, 1.0, 1.0)
private let negOne = nqScale(-1.0, one)
private let half = nqScale(0.5, one)
private let negHalf = nqScale(-1.0, half)
private let two = nqScale(2.0, one)
private let negTwo = nqScale(-1.0, two)
private let three = nqScale(3.0, one)
private let negThree = nqScale(-1.0, three)
private let four = nqScale(4.0, one)
private let six = nqScale(6.0, one)

private class DualTest {
    let name: String
    let x: Array<NqNumber>
    let fnDual: (NdqNumber) -> NdqNumber
    let fn: (NqNumber) -> NqNumber
    let dFn: (NqNumber) -> NqNumber
    public init(name: String, x: Array<NqNumber>, fnDual: (NdqNumber) -> NdqNumber, fn: (NqNumber) -> NqNumber,
        dFn: (NqNumber) -> NqNumber) {
        this.name = name
        this.x = x
        this.fnDual = fnDual
        this.fn = fn
        this.dFn = dFn
    }
}

private let dualTests: Array<DualTest> = [
    DualTest(
        "exp",
        [nqNaN(), nqInf(), negThree, negTwo, negOne, negHalf, negZeroQuat, NqNumber(0.0, 0.0, 0.0, 0.0), half, one, two,
            three],
        ndqExp,
        nqExp,
        dExp,
    ),
    DualTest(
        "log",
        [nqNaN(), nqInf(), negThree, negTwo, negOne, negHalf, negZeroQuat, NqNumber(0.0, 0.0, 0.0, 0.0), half, one, two,
            three],
        ndqLog,
        nqLog,
        dLog,
    ),
    DualTest(
        "inv",
        [nqNaN(), nqInf(), negThree, negTwo, negOne, negHalf, negZeroQuat, NqNumber(0.0, 0.0, 0.0, 0.0), half, one, two,
            three],
        ndqInv,
        nqInv,
        dInv,
    ),
    DualTest(
        "sqrt",
        [nqNaN(), nqInf(), negThree, negTwo, negOne, negHalf, negZeroQuat, NqNumber(0.0, 0.0, 0.0, 0.0), half, one, two,
            three],
        ndqSqrt,
        nqSqrt,
        dSqrt,
    )
]

@Test
public func testNumber(): Unit {
    let tol = 1e-15
    for (test in dualTests) {
        for (x in test.x) {
            let fxDual = test.fnDual(NdqNumber(x, oneReal))
            let fx = test.fn(x)
            let dFx = test.dFn(x)
            @Expect(same(fxDual.real, fx, tol))
            @Expect(same(fxDual.dual, dFx, tol))
        }
    }
}

private let invTests: Array<NdqNumber> = [
    NdqNumber(NqNumber(1.0, 0.0, 0.0, 0.0), NqNumber(0.0, 0.0, 0.0, 0.0)),
    NdqNumber(
        NqNumber(1.0, 0.0, 0.0, 0.0),
        NqNumber(1.0, 0.0, 0.0, 0.0)
    ),
    NdqNumber(
        NqNumber(0.0, 1.0, 0.0, 0.0),
        NqNumber(1.0, 0.0, 0.0, 0.0)
    ),
    NdqNumber(
        NqNumber(1.0, 0.0, 0.0, 0.0),
        NqNumber(1.0, 1.0, 0.0, 0.0)
    ),
    NdqNumber(
        NqNumber(1.0, 1.0, 0.0, 0.0),
        NqNumber(1.0, 1.0, 0.0, 0.0)
    ),
    NdqNumber(
        NqNumber(1.0, 10.0, 0.0, 0.0),
        NqNumber(1.0, 5.0, 0.0, 0.0)
    ),
    NdqNumber(
        NqNumber(10.0, 1.0, 0.0, 0.0),
        NqNumber(5.0, 1.0, 0.0, 0.0)
    ),
    NdqNumber(
        NqNumber(1.0, 0.0, 0.0, 0.0),
        NqNumber(1.0, 1.0, 0.0, 1.0)
    ),
    NdqNumber(
        NqNumber(12.0, 1.0, 0.0, 0.0),
        NqNumber(1.0, 1.0, 0.0, 0.0)
    ),
    NdqNumber(
        NqNumber(12.0, 1.0, 3.0, 0.0),
        NqNumber(0.0, 0.0, 0.0, 0.0)
    )
]

@Test
public func testInv(): Unit {
    let tol = 1e-15
    for (x in invTests) {
        let got = ndqMul(x, ndqInv(x))
        let want = NdqNumber(NqNumber(1.0, 0.0, 0.0, 0.0), NqNumber(0.0, 0.0, 0.0, 0.0))
        @Expect(sameDual(got, want, tol))
    }
}

private class PowRealTest {
    let d: NdqNumber
    let p: Float64
    let want: NdqNumber
    public init(d: NdqNumber, p: Float64, want: NdqNumber) {
        this.d = d
        this.p = p
        this.want = want
    }
}

private let powRealTests: Array<PowRealTest> = [
    // PowReal(NaN+xϵ, ±0) = 1+NaNϵ for any x
    PowRealTest(NdqNumber(nqNaN(), zeroQuat), 0.0, NdqNumber(oneReal, nqNaN())),
    PowRealTest(NdqNumber(nqNaN(), zeroQuat), negZero, NdqNumber(oneReal, nqNaN())),
    PowRealTest(NdqNumber(nqNaN(), one), 0.0, NdqNumber(oneReal, nqNaN())),
    PowRealTest(NdqNumber(nqNaN(), two), negZero, NdqNumber(oneReal, nqNaN())),
    PowRealTest(NdqNumber(nqNaN(), three), 0.0, NdqNumber(oneReal, nqNaN())),
    PowRealTest(NdqNumber(nqNaN(), one), negZero, NdqNumber(oneReal, nqNaN())),
    PowRealTest(NdqNumber(nqNaN(), two), 0.0, NdqNumber(oneReal, nqNaN())),
    PowRealTest(NdqNumber(nqNaN(), three), negZero, NdqNumber(oneReal, nqNaN())),

    // PowReal(x, ±0) = 1 for any x
    PowRealTest(NdqNumber(zeroQuat, zeroQuat), 0.0, NdqNumber(oneReal, zeroQuat)),
    PowRealTest(NdqNumber(negZeroQuat, zeroQuat), negZero, NdqNumber(oneReal, zeroQuat)),
    PowRealTest(NdqNumber(nqInf(), zeroQuat), 0.0, NdqNumber(oneReal, zeroQuat)),
    PowRealTest(NdqNumber(nqInf(), zeroQuat), negZero, NdqNumber(oneReal, zeroQuat)),
    PowRealTest(NdqNumber(zeroQuat, one), 0.0, NdqNumber(oneReal, zeroQuat)),
    PowRealTest(NdqNumber(negZeroQuat, one), negZero, NdqNumber(oneReal, zeroQuat)),
    PowRealTest(NdqNumber(nqInf(), one), 0.0, NdqNumber(oneReal, zeroQuat)),
    PowRealTest(NdqNumber(nqInf(), one), negZero, NdqNumber(oneReal, zeroQuat)),

    // PowReal(1+xϵ, y) = (1+xyϵ) for any y
    PowRealTest(NdqNumber(oneReal, zeroQuat), 0.0, NdqNumber(oneReal, zeroQuat)),
    PowRealTest(NdqNumber(oneReal, zeroQuat), 1.0, NdqNumber(oneReal, zeroQuat)),
    PowRealTest(NdqNumber(oneReal, zeroQuat), 2.0, NdqNumber(oneReal, zeroQuat)),
    PowRealTest(NdqNumber(oneReal, zeroQuat), 3.0, NdqNumber(oneReal, zeroQuat)),
    PowRealTest(NdqNumber(oneReal, one), 0.0, NdqNumber(oneReal, zeroQuat)),
    PowRealTest(NdqNumber(oneReal, one), 1.0, NdqNumber(oneReal, one)),
    PowRealTest(NdqNumber(oneReal, one), 2.0, NdqNumber(oneReal, two)),
    PowRealTest(NdqNumber(oneReal, one), 3.0, NdqNumber(oneReal, three)),
    PowRealTest(NdqNumber(oneReal, two), 0.0, NdqNumber(oneReal, zeroQuat)),
    PowRealTest(NdqNumber(oneReal, two), 1.0, NdqNumber(oneReal, two)),
    PowRealTest(NdqNumber(oneReal, two), 2.0, NdqNumber(oneReal, four)),
    PowRealTest(NdqNumber(oneReal, two), 3.0, NdqNumber(oneReal, six)),

    // PowReal(x, 1) = x for any x
    PowRealTest(NdqNumber(zeroQuat, zeroQuat), 1.0, NdqNumber(zeroQuat, zeroQuat)),
    PowRealTest(NdqNumber(negZeroQuat, zeroQuat), 1.0, NdqNumber(negZeroQuat, zeroQuat)),
    PowRealTest(NdqNumber(zeroQuat, one), 1.0, NdqNumber(zeroQuat, one)),
    PowRealTest(NdqNumber(negZeroQuat, one), 1.0, NdqNumber(negZeroQuat, one)),
    PowRealTest(NdqNumber(nqNaN(), zeroQuat), 1.0, NdqNumber(nqNaN(), zeroQuat)),
    PowRealTest(NdqNumber(nqNaN(), one), 1.0, NdqNumber(nqNaN(), one)),
    PowRealTest(NdqNumber(nqNaN(), two), 1.0, NdqNumber(nqNaN(), two)),

    // PowReal(NaN+xϵ, y) = NaN+NaNϵ
    PowRealTest(NdqNumber(nqNaN(), zeroQuat), 2.0, NdqNumber(nqNaN(), nqNaN())),
    PowRealTest(NdqNumber(nqNaN(), zeroQuat), 3.0, NdqNumber(nqNaN(), nqNaN())),
    PowRealTest(NdqNumber(nqNaN(), one), 2.0, NdqNumber(nqNaN(), nqNaN())),
    PowRealTest(NdqNumber(nqNaN(), one), 3.0, NdqNumber(nqNaN(), nqNaN())),
    PowRealTest(NdqNumber(nqNaN(), two), 2.0, NdqNumber(nqNaN(), nqNaN())),
    PowRealTest(NdqNumber(nqNaN(), two), 3.0, NdqNumber(nqNaN(), nqNaN())),

    // PowReal(x, NaN) = NaN+NaNϵ
    PowRealTest(NdqNumber(zeroQuat, zeroQuat), Float64.NaN, NdqNumber(nqNaN(), nqNaN())),
    PowRealTest(NdqNumber(two, zeroQuat), Float64.NaN, NdqNumber(nqNaN(), nqNaN())),
    PowRealTest(NdqNumber(three, zeroQuat), Float64.NaN, NdqNumber(nqNaN(), nqNaN())),
    PowRealTest(NdqNumber(zeroQuat, one), Float64.NaN, NdqNumber(nqNaN(), nqNaN())),
    PowRealTest(NdqNumber(two, one), Float64.NaN, NdqNumber(nqNaN(), nqNaN())),
    PowRealTest(NdqNumber(three, one), Float64.NaN, NdqNumber(nqNaN(), nqNaN())),
    PowRealTest(NdqNumber(zeroQuat, two), Float64.NaN, NdqNumber(nqNaN(), nqNaN())),
    PowRealTest(NdqNumber(two, two), Float64.NaN, NdqNumber(nqNaN(), nqNaN())),
    PowRealTest(NdqNumber(three, two), Float64.NaN, NdqNumber(nqNaN(), nqNaN())),

    // Handled by quat.Pow tests:
    //
    // Pow(±0, y) = ±Inf for y an odd integer < 0
    // Pow(±0, -Inf) = +Inf
    // Pow(±0, +Inf) = +0
    // Pow(±0, y) = +Inf for finite y < 0 and not an odd integer
    // Pow(±0, y) = ±0 for y an odd integer > 0
    // Pow(±0, y) = +0 for finite y > 0 and not an odd integer
    // Pow(-1, ±Inf) = 1

    // PowReal(x+0ϵ, +Inf) = +Inf+NaNϵ for |x| > 1
    PowRealTest(NdqNumber(two, zeroQuat), Float64.Inf, NdqNumber(nqInf(), nqNaN())),
    PowRealTest(NdqNumber(three, zeroQuat), Float64.Inf, NdqNumber(nqInf(), nqNaN())),

    // PowReal(x+yϵ, +Inf) = +Inf for |x| > 1
    PowRealTest(NdqNumber(two, one), Float64.Inf, NdqNumber(nqInf(), nqInf())),
    PowRealTest(NdqNumber(three, one), Float64.Inf, NdqNumber(nqInf(), nqInf())),
    PowRealTest(NdqNumber(two, two), Float64.Inf, NdqNumber(nqInf(), nqInf())),
    PowRealTest(NdqNumber(three, two), Float64.Inf, NdqNumber(nqInf(), nqInf())),

    // PowReal(x, -Inf) = +0+NaNϵ for |x| > 1
    PowRealTest(NdqNumber(two, zeroQuat), -Float64.Inf, NdqNumber(zeroQuat, nqNaN())),
    PowRealTest(NdqNumber(three, zeroQuat), -Float64.Inf, NdqNumber(zeroQuat, nqNaN())),
    PowRealTest(NdqNumber(two, one), -Float64.Inf, NdqNumber(zeroQuat, nqNaN())),
    PowRealTest(NdqNumber(three, one), -Float64.Inf, NdqNumber(zeroQuat, nqNaN())),
    PowRealTest(NdqNumber(two, two), -Float64.Inf, NdqNumber(zeroQuat, nqNaN())),
    PowRealTest(NdqNumber(three, two), -Float64.Inf, NdqNumber(zeroQuat, nqNaN())),

    // PowReal(x+yϵ, +Inf) = +0+NaNϵ for |x| < 1
    PowRealTest(NdqNumber(nqScale(0.1, one), zeroQuat), Float64.Inf, NdqNumber(zeroQuat, nqNaN())),
    PowRealTest(NdqNumber(nqScale(0.1, one), nqScale(0.1, one)), Float64.Inf, NdqNumber(zeroQuat, nqNaN())),
    PowRealTest(NdqNumber(nqScale(0.2, one), nqScale(0.2, one)), Float64.Inf, NdqNumber(zeroQuat, nqNaN())),
    PowRealTest(NdqNumber(nqScale(0.5, one), nqScale(0.5, one)), Float64.Inf, NdqNumber(zeroQuat, nqNaN())),

    // PowReal(x+0ϵ, -Inf) = +Inf+NaNϵ for |x| < 1
    PowRealTest(NdqNumber(nqScale(0.1, one), zeroQuat), -Float64.Inf, NdqNumber(nqInf(), nqNaN())),
    PowRealTest(NdqNumber(nqScale(0.2, one), zeroQuat), -Float64.Inf, NdqNumber(nqInf(), nqNaN())),

    // PowReal(x, -Inf) = +Inf-Infϵ for |x| < 1
    PowRealTest(NdqNumber(nqScale(0.1, one), nqScale(0.1, one)), -Float64.Inf, NdqNumber(nqInf(), nqInf())),
    PowRealTest(NdqNumber(nqScale(0.2, one), nqScale(0.1, one)), -Float64.Inf, NdqNumber(nqInf(), nqInf())),
    PowRealTest(NdqNumber(nqScale(0.1, one), nqScale(0.2, one)), -Float64.Inf, NdqNumber(nqInf(), nqInf())),
    PowRealTest(NdqNumber(nqScale(0.2, one), nqScale(0.2, one)), -Float64.Inf, NdqNumber(nqInf(), nqInf())),
    PowRealTest(NdqNumber(nqScale(0.1, one), one), -Float64.Inf, NdqNumber(nqInf(), nqInf())),
    PowRealTest(NdqNumber(nqScale(0.2, one), one), -Float64.Inf, NdqNumber(nqInf(), nqInf())),
    PowRealTest(NdqNumber(nqScale(0.1, one), two), -Float64.Inf, NdqNumber(nqInf(), nqInf())),
    PowRealTest(NdqNumber(nqScale(0.2, one), two), -Float64.Inf, NdqNumber(nqInf(), nqInf()))

    // Handled by quat.Pow tests:
    //
    // Pow(+Inf, y) = +Inf for y > 0
    // Pow(+Inf, y) = +0 for y < 0
    // Pow(-Inf, y) = Pow(-0, -y)
]

@Test
public func testPowReal(): Unit {
    let tol = 1e-15
    for (test in powRealTests) {
        let got = ndqPowReal(test.d, test.p)
        @Expect(sameDual(got, test.want, tol))
    }
}

private func sameDual(a: NdqNumber, b: NdqNumber, tol: Float64): Bool {
    return same(a.real, b.real, tol) && same(a.dual, b.dual, tol)
}

private func same(a: NqNumber, b: NqNumber, tol: Float64): Bool {
    return (nqIsNaN(a) && nqIsNaN(b)) || (nqIsInf(a) && nqIsInf(b)) || equalApprox(a, b, tol)
}

private func equalApprox(a: NqNumber, b: NqNumber, tol: Float64): Bool {
    return equalWithinAbsOrRel(a.real, b.real, tol, tol) && equalWithinAbsOrRel(a.imag, b.imag, tol, tol) &&
        equalWithinAbsOrRel(a.jmag, b.jmag, tol, tol) && equalWithinAbsOrRel(a.kmag, b.kmag, tol, tol)
}
