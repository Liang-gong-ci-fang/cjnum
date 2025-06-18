package cjnum.num.hyperdual

import cjnum.floats.scalar.*
import cjnum.util.*
import std.math.*
import std.unittest.*
import std.unittest.testmacro.*

// TODO var formatTests = []struct {

// TODO func TestFormat(t *testing.T) {

// First derivatives:

private func dSin(x: Float64): Float64 {
    return cos(x)
}

private func dCos(x: Float64): Float64 {
    return -sin(x)
}

private func dTan(x: Float64): Float64 {
    return sec(x) * sec(x)
}

private func dAsin(x: Float64): Float64 {
    return 1.0 / nanSqrt(1.0 - x * x)
}

private func dAcos(x: Float64): Float64 {
    return -1.0 / nanSqrt(1.0 - x * x)
}

private func dAtan(x: Float64): Float64 {
    return 1.0 / (1.0 + x * x)
}

private func dSinh(x: Float64): Float64 {
    return cosh(x)
}

private func dCosh(x: Float64): Float64 {
    return sinh(x)
}

private func dTanh(x: Float64): Float64 {
    return sech(x) * sech(x)
}

private func dAsinh(x: Float64): Float64 {
    return 1.0 / sqrt(x * x + 1.0)
}

private func dAcosh(x: Float64): Float64 {
    return 1.0 / (nanSqrt(x - 1.0) * nanSqrt(x + 1.0))
}

private func dAtanh(x: Float64): Float64 {
    if (abs(x) == 1.0) {
        return Float64.NaN
    } else if (x.isInf()) {
        return negZero
    }
    return 1.0 / (1.0 - x * x)
}

private func dExp(x: Float64): Float64 {
    return exp(x)
}

private func dLog(x: Float64): Float64 {
    if (x < 0.0) {
        return Float64.NaN
    }
    return 1.0 / x
}

private func dSqrt(x: Float64): Float64 {
    // For whatever reason, nanSqrt(-0) returns -0.
    // In this case, that is clearly a wrong approach.
    if (x < 0.0) {
        return Float64.NaN
    } else if (x == 0.0) {
        return Float64.Inf
    }
    return 0.5 / sqrt(x)
}

private func dInv(x: Float64): Float64 {
    return -1.0 / (x * x)
}

// Second derivatives:

private func d2Sin(x: Float64): Float64 {
    return -sin(x)
}

private func d2Cos(x: Float64): Float64 {
    return -cos(x)
}

private func d2Tan(x: Float64): Float64 {
    return 2.0 * tan(x) * sec(x) * sec(x)
}

private func d2Asin(x: Float64): Float64 {
    return x / pow(1.0 - x * x, 1.5)
}

private func d2Acos(x: Float64): Float64 {
    return -x / pow(1.0 - x * x, 1.5)
}

private func d2Atan(x: Float64): Float64 {
    return -2.0 * x / ((x * x + 1.0) * (x * x + 1.0))
}

private func d2Sinh(x: Float64): Float64 {
    return nanSinh(x)
}

private func d2Cosh(x: Float64): Float64 {
    return nanCosh(x)
}

private func d2Tanh(x: Float64): Float64 {
    return -2.0 * nanTanh(x) * sech(x) * sech(x)
}

private func d2Asinh(x: Float64): Float64 {
    return -x / pow((x * x + 1.0), 1.5)
}

private func d2Acosh(x: Float64): Float64 {
    return -x / (pow(x - 1.0, 1.5) * pow(x + 1.0, 1.5))
}

private func d2Atanh(x: Float64): Float64 {
    return 2.0 * x / ((1.0 - x * x) * (1.0 - x * x))
}

private func d2Exp(x: Float64): Float64 {
    return exp(x)
}

private func d2Log(x: Float64): Float64 {
    if (x < 0.0) {
        return Float64.NaN
    }
    return -1.0 / (x * x)
}

private func d2Sqrt(x: Float64): Float64 {
    if (x == Float64.Inf) {
        return 0.0
    } else if (x == -Float64.Inf) {
        return Float64.NaN
    }
    return -0.25 * pow(x, -1.5)
}

func d2Inv(x: Float64): Float64 {
    return 2.0 / (x * x * x)
}

// Helpers:

private func sec(x: Float64): Float64 {
    return 1.0 / cos(x)
}

private func sech(x: Float64): Float64 {
    return 1.0 / cosh(x)
}

private class HyperdualTest {
    let name: String
    let x: Array<Float64>
    let fnHyperdual: (x: NhNumber) -> NhNumber
    let fn: (x: Float64) -> Float64
    let dFn: (x: Float64) -> Float64
    let d2Fn: (x: Float64) -> Float64
    public init(name: String, x: Array<Float64>, fnHyperdual: (x: NhNumber) -> NhNumber, fn: (x: Float64) -> Float64,
        dFn: (x: Float64) -> Float64, d2Fn: (x: Float64) -> Float64) {
        this.name = name
        this.x = x
        this.fnHyperdual = fnHyperdual
        this.fn = fn
        this.dFn = dFn
        this.d2Fn = d2Fn
    }
}

private let hyperdualTests: Array<HyperdualTest> = [
    HyperdualTest(
        "sin",
        [Float64.NaN, -Float64.Inf, -3.0, -2.0, -1.0, -0.5, negZero, 0.0, 0.5, 1.0, 2.0, 3.0, Float64.Inf],
        nhSin,
        sin,
        dSin,
        d2Sin,
    ),
    HyperdualTest(
        "cos",
        [Float64.NaN, -Float64.Inf, -3.0, -2.0, -1.0, -0.5, negZero, 0.0, 0.5, 1.0, 2.0, 3.0, Float64.Inf],
        nhCos,
        cos,
        dCos,
        d2Cos,
    ),
    HyperdualTest(
        "tan",
        [Float64.NaN, -Float64.Inf, -3.0, -2.0, -1.0, -0.5, negZero, 0.0, 0.5, 1.0, 2.0, 3.0, Float64.Inf],
        nhTan,
        tan,
        dTan,
        d2Tan,
    ),
    HyperdualTest(
        "sinh",
        [Float64.NaN, -Float64.Inf, -3.0, -2.0, -1.0, -0.5, negZero, 0.0, 0.5, 1.0, 2.0, 3.0, Float64.Inf],
        nhSinh,
        sinh,
        dSinh,
        d2Sinh,
    ),
    HyperdualTest(
        "cosh",
        [Float64.NaN, -Float64.Inf, -3.0, -2.0, -1.0, -0.5, negZero, 0.0, 0.5, 1.0, 2.0, 3.0, Float64.Inf],
        nhCosh,
        cosh,
        dCosh,
        d2Cosh,
    ),
    HyperdualTest(
        "tanh",
        [Float64.NaN, -Float64.Inf, -3.0, -2.0, -1.0, -0.5, negZero, 0.0, 0.5, 1.0, 2.0, 3.0, Float64.Inf],
        nhTanh,
        nanTanh,
        dTanh,
        d2Tanh,
    ),
    HyperdualTest(
        "asin",
        [Float64.NaN, -Float64.Inf, -3.0, -2.0, -1.0, -0.5, negZero, 0.0, 0.5, 1.0, 2.0, 3.0, Float64.Inf],
        nhAsin,
        nanAsin,
        dAsin,
        d2Asin,
    ),
    HyperdualTest(
        "acos",
        [Float64.NaN, -Float64.Inf, -3.0, -2.0, -1.0, -0.5, negZero, 0.0, 0.5, 1.0, 2.0, 3.0, Float64.Inf],
        nhAcos,
        nanAcos,
        dAcos,
        d2Acos,
    ),
    HyperdualTest(
        "atan",
        [Float64.NaN, -Float64.Inf, -3.0, -2.0, -1.0, -0.5, negZero, 0.0, 0.5, 1.0, 2.0, 3.0, Float64.Inf],
        nhAtan,
        atan,
        dAtan,
        d2Atan,
    ),
    HyperdualTest(
        "asinh",
        [Float64.NaN, -Float64.Inf, -3.0, -2.0, -1.0, -0.5, negZero, 0.0, 0.5, 1.0, 2.0, 3.0, Float64.Inf],
        nhAsinh,
        asinh,
        dAsinh,
        d2Asinh,
    ),
    HyperdualTest(
        "acosh",
        [Float64.NaN, -Float64.Inf, -3.0, -2.0, -1.0, -0.5, negZero, 0.0, 0.5, 1.0, 2.0, 3.0, Float64.Inf],
        nhAcosh,
        nanAcosh,
        dAcosh,
        d2Acosh,
    ),
    HyperdualTest(
        "atanh",
        [Float64.NaN, -Float64.Inf, -3.0, -2.0, -1.0, -0.5, negZero, 0.0, 0.5, 1.0, 2.0, 3.0, Float64.Inf],
        nhAtanh,
        nanAtanh,
        dAtanh,
        d2Atanh,
    ),
    HyperdualTest(
        "exp",
        [Float64.NaN, -Float64.Inf, -3.0, -2.0, -1.0, -0.5, negZero, 0.0, 0.5, 1.0, 2.0, 3.0, Float64.Inf],
        nhExp,
        exp,
        dExp,
        d2Exp,
    ),
    HyperdualTest(
        "log",
        [Float64.NaN, -Float64.Inf, -3.0, -2.0, -1.0, -0.5, negZero, 0.0, 0.5, 1.0, 2.0, 3.0, Float64.Inf],
        nhLog,
        log,
        dLog,
        d2Log,
    ),
    HyperdualTest(
        "inv",
        [Float64.NaN, -Float64.Inf, -3.0, -2.0, -1.0, -0.5, negZero, 0.0, 0.5, 1.0, 2.0, 3.0, Float64.Inf],
        nhInv,
        {x: Float64 => 1.0 / x},
        dInv,
        d2Inv,
    ),
    HyperdualTest(
        "sqrt",
        [Float64.NaN, -Float64.Inf, -3.0, -2.0, -1.0, -0.5, negZero, 0.0, 0.5, 1.0, 2.0, 3.0, Float64.Inf],
        nhSqrt,
        nanSqrt,
        dSqrt,
        d2Sqrt,
    ),
    HyperdualTest(
        "Fike example fn",
        [1.0, 2.0, 3.0, 4.0, 5.0],
        {
            x: NhNumber => nhMul(
                nhExp(x),
                nhInv(
                    nhSqrt(
                        nhAdd(
                            nhPowReal(nhSin(x), 3.0),
                            nhPowReal(nhCos(x), 3.0)
                        )
                    )
                )
            )
        },
        {x: Float64 => exp(x) / nanSqrt(pow(sin(x), 3.0) + pow(cos(x), 3.0))},
        {
            x: Float64 => exp(x) * (3.0 * cos(x) + 5.0 * cos(3.0 * x) + 9.0 * sin(x) + sin(3.0 * x)) / (8.0 * pow(
                pow(sin(x), 3.0) + pow(cos(x), 3.0), 1.5))
        },
        {
            x: Float64 => exp(x) * (130.0 - 12.0 * cos(2.0 * x) + 30.0 * cos(4.0 * x) + 12.0 * cos(6.0 * x) - 111.0 *
                sin(2.0 * x) + 48.0 * sin(4.0 * x) + 5.0 * sin(6.0 * x)) / (64.0 * pow(
                pow(sin(x), 3.0) + pow(cos(x), 3.0),
                2.5
            ))
        },
    )
]

@Test
public func testHyperdual(): Unit {
    let tol = 1e-14
    for (test in hyperdualTests) {
        for (x in test.x) {
            let fxHyperdual = test.fnHyperdual(NhNumber(x, 1.0, 1.0, 0.0))
            let fx = test.fn(x)
            let dFx = test.dFn(x)
            let d2Fx = test.d2Fn(x)
            @Expect(same(fxHyperdual.real, fx, tol))
            @Expect(same(fxHyperdual.e1mag, dFx, tol))
            @Expect(same(fxHyperdual.e1mag, fxHyperdual.e2mag, tol))
            @Expect(same(fxHyperdual.e1E2mag, d2Fx, tol))
        }
    }
}

private class PowRealTest {
    let d: NhNumber
    let p: Float64
    let want: NhNumber
    public init(d: NhNumber, p: Float64, want: NhNumber) {
        this.d = d
        this.p = p
        this.want = want
    }
}

private let powRealTests = [
    // PowReal(NaN+xϵ₁+yϵ₂, ±0) = 1+NaNϵ₁+NaNϵ₂+NaNϵ₁ϵ₂ for any x and y
    PowRealTest(NhNumber(Float64.NaN, 0.0, 0.0, 0.0), 0.0, NhNumber(1.0, Float64.NaN, Float64.NaN, Float64.NaN)),
    // PowReal(x, NaN) = NaN+NaNϵ
    PowRealTest(NhNumber(Float64.NaN, 0.0, 0.0, 0.0), negZero, NhNumber(1.0, Float64.NaN, Float64.NaN, Float64.NaN)),
    PowRealTest(NhNumber(Float64.NaN, 1.0, 1.0, 0.0), 0.0, NhNumber(1.0, Float64.NaN, Float64.NaN, Float64.NaN)),
    PowRealTest(NhNumber(Float64.NaN, 2.0, 2.0, 0.0), negZero, NhNumber(1.0, Float64.NaN, Float64.NaN, Float64.NaN)),
    PowRealTest(NhNumber(Float64.NaN, 3.0, 3.0, 0.0), 0.0, NhNumber(1.0, Float64.NaN, Float64.NaN, Float64.NaN)),
    PowRealTest(NhNumber(Float64.NaN, 1.0, 1.0, 0.0), negZero, NhNumber(1.0, Float64.NaN, Float64.NaN, Float64.NaN)),
    PowRealTest(NhNumber(Float64.NaN, 2.0, 2.0, 0.0), 0.0, NhNumber(1.0, Float64.NaN, Float64.NaN, Float64.NaN)),
    PowRealTest(NhNumber(Float64.NaN, 3.0, 3.0, 0.0), negZero, NhNumber(1.0, Float64.NaN, Float64.NaN, Float64.NaN)),
    PowRealTest(NhNumber(Float64.NaN, 2.0, 3.0, 0.0), 0.0, NhNumber(1.0, Float64.NaN, Float64.NaN, Float64.NaN)),
    PowRealTest(NhNumber(Float64.NaN, 2.0, 3.0, 0.0), negZero, NhNumber(1.0, Float64.NaN, Float64.NaN, Float64.NaN)),

    // PowReal(x, ±0) = 1 for any x
    PowRealTest(NhNumber(0.0, 0.0, 0.0, 0.0), 0.0, NhNumber(1.0, 0.0, 0.0, 0.0)),
    PowRealTest(NhNumber(Float64.Inf, 0.0, 0.0, 0.0), 0.0, NhNumber(1.0, 0.0, 0.0, 0.0)),
    PowRealTest(NhNumber(-Float64.Inf, 0.0, 0.0, 0.0), negZero, NhNumber(1.0, 0.0, 0.0, 0.0)),
    PowRealTest(NhNumber(0.0, 1.0, 1.0, 0.0), 0.0, NhNumber(1.0, 0.0, 0.0, 0.0)),
    PowRealTest(NhNumber(Float64.Inf, 1.0, 1.0, 0.0), 0.0, NhNumber(1.0, 0.0, 0.0, 0.0)),
    PowRealTest(NhNumber(-Float64.Inf, 1.0, 1.0, 0.0), negZero, NhNumber(1.0, 0.0, 0.0, 0.0)),
    // These two satisfy the claim above, but the sign of zero is negative. Do we care?
    PowRealTest(NhNumber(negZero, 0.0, 0.0, 0.0), negZero, NhNumber(1.0, negZero, negZero, 0.0)),
    PowRealTest(NhNumber(negZero, 1.0, 1.0, 0.0), negZero, NhNumber(1.0, negZero, negZero, 0.0)),

    // PowReal(1+xϵ₁+yϵ₂, z) = 1+xzϵ₁+yzϵ₂+2xyzϵ₁ϵ₂ for any z
    PowRealTest(NhNumber(1.0, 0.0, 0.0, 0.0), 0.0, NhNumber(1.0, 0.0, 0.0, 0.0)),
    PowRealTest(NhNumber(1.0, 0.0, 0.0, 0.0), 1.0, NhNumber(1.0, 0.0, 0.0, 0.0)),
    PowRealTest(NhNumber(1.0, 0.0, 0.0, 0.0), 2.0, NhNumber(1.0, 0.0, 0.0, 0.0)),
    PowRealTest(NhNumber(1.0, 0.0, 0.0, 0.0), 3.0, NhNumber(1.0, 0.0, 0.0, 0.0)),
    PowRealTest(NhNumber(1.0, 1.0, 1.0, 0.0), 0.0, NhNumber(1.0, 0.0, 0.0, 0.0)),
    PowRealTest(NhNumber(1.0, 1.0, 1.0, 0.0), 1.0, NhNumber(1.0, 1.0, 1.0, 0.0)),
    PowRealTest(NhNumber(1.0, 1.0, 1.0, 0.0), 2.0, NhNumber(1.0, 2.0, 2.0, 2.0)),
    PowRealTest(NhNumber(1.0, 1.0, 1.0, 0.0), 3.0, NhNumber(1.0, 3.0, 3.0, 6.0)),
    PowRealTest(NhNumber(1.0, 2.0, 2.0, 0.0), 0.0, NhNumber(1.0, 0.0, 0.0, 0.0)),
    PowRealTest(NhNumber(1.0, 2.0, 2.0, 0.0), 1.0, NhNumber(1.0, 2.0, 2.0, 0.0)),
    PowRealTest(NhNumber(1.0, 2.0, 2.0, 0.0), 2.0, NhNumber(1.0, 4.0, 4.0, 8.0)),
    PowRealTest(NhNumber(1.0, 2.0, 2.0, 0.0), 3.0, NhNumber(1.0, 6.0, 6.0, 24.0)),
    PowRealTest(NhNumber(1.0, 1.0, 2.0, 0.0), 0.0, NhNumber(1.0, 0.0, 0.0, 0.0)),
    PowRealTest(NhNumber(1.0, 1.0, 2.0, 0.0), 1.0, NhNumber(1.0, 1.0, 2.0, 0.0)),
    PowRealTest(NhNumber(1.0, 1.0, 2.0, 0.0), 2.0, NhNumber(1.0, 2.0, 4.0, 4.0)),
    PowRealTest(NhNumber(1.0, 1.0, 2.0, 0.0), 3.0, NhNumber(1.0, 3.0, 6.0, 12.0)),

    // PowReal(NaN+xϵ₁+yϵ₂, 1) = NaN+xϵ₁+yϵ₂+NaNϵ₁ϵ₂ for any x
    PowRealTest(NhNumber(Float64.NaN, 0.0, 0.0, 0.0), 1.0, NhNumber(Float64.NaN, 0.0, 0.0, Float64.NaN)),
    PowRealTest(NhNumber(Float64.NaN, 1.0, 1.0, 0.0), 1.0, NhNumber(Float64.NaN, 1.0, 1.0, Float64.NaN)),
    PowRealTest(NhNumber(Float64.NaN, 2.0, 2.0, 0.0), 1.0, NhNumber(Float64.NaN, 2.0, 2.0, Float64.NaN)),
    PowRealTest(NhNumber(Float64.NaN, 1.0, 2.0, 0.0), 1.0, NhNumber(Float64.NaN, 1.0, 2.0, Float64.NaN)),

    // PowReal(x, 1) = x for any x
    PowRealTest(NhNumber(0.0, 0.0, 0.0, 0.0), 1.0, NhNumber(0.0, 0.0, 0.0, 0.0)),
    PowRealTest(NhNumber(negZero, 0.0, 0.0, 0.0), 1.0, NhNumber(negZero, 0.0, 0.0, 0.0)),
    PowRealTest(NhNumber(0.0, 1.0, 1.0, 0.0), 1.0, NhNumber(0.0, 1.0, 1.0, 0.0)),
    PowRealTest(NhNumber(negZero, 1.0, 1.0, 0.0), 1.0, NhNumber(negZero, 1.0, 1.0, 0.0)),
    PowRealTest(NhNumber(0.0, 1.0, 2.0, 0.0), 1.0, NhNumber(0.0, 1.0, 2.0, 0.0)),
    PowRealTest(NhNumber(negZero, 1.0, 2.0, 0.0), 1.0, NhNumber(negZero, 1.0, 2.0, 0.0)),

    // PowReal(NaN+xϵ₁+xϵ₂, y) = NaN+NaNϵ₁+NaNϵ₂+NaNϵ₁ϵ₂
    PowRealTest(NhNumber(Float64.NaN, 0.0, 0.0, 0.0), 2.0, NhNumber(Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN)),
    PowRealTest(NhNumber(Float64.NaN, 0.0, 0.0, 0.0), 3.0, NhNumber(Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN)),
    PowRealTest(NhNumber(Float64.NaN, 1.0, 1.0, 0.0), 2.0, NhNumber(Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN)),
    PowRealTest(NhNumber(Float64.NaN, 1.0, 1.0, 0.0), 3.0, NhNumber(Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN)),
    PowRealTest(NhNumber(Float64.NaN, 2.0, 2.0, 0.0), 2.0, NhNumber(Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN)),
    PowRealTest(NhNumber(Float64.NaN, 2.0, 2.0, 0.0), 3.0, NhNumber(Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN)),
    PowRealTest(NhNumber(Float64.NaN, 1.0, 2.0, 0.0), 2.0, NhNumber(Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN)),
    PowRealTest(NhNumber(Float64.NaN, 1.0, 2.0, 0.0), 3.0, NhNumber(Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN)),

    // PowReal(x, NaN) = NaN+NaNϵ₁+NaNϵ₂+NaNϵ₁ϵ₂
    PowRealTest(NhNumber(0.0, 0.0, 0.0, 0.0), Float64.NaN, NhNumber(Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN)),
    PowRealTest(NhNumber(2.0, 0.0, 0.0, 0.0), Float64.NaN, NhNumber(Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN)),
    PowRealTest(NhNumber(3.0, 0.0, 0.0, 0.0), Float64.NaN, NhNumber(Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN)),
    PowRealTest(NhNumber(0.0, 1.0, 1.0, 0.0), Float64.NaN, NhNumber(Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN)),
    PowRealTest(NhNumber(2.0, 1.0, 1.0, 0.0), Float64.NaN, NhNumber(Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN)),
    PowRealTest(NhNumber(3.0, 1.0, 1.0, 0.0), Float64.NaN, NhNumber(Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN)),
    PowRealTest(NhNumber(0.0, 2.0, 2.0, 0.0), Float64.NaN, NhNumber(Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN)),
    PowRealTest(NhNumber(2.0, 2.0, 2.0, 0.0), Float64.NaN, NhNumber(Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN)),
    PowRealTest(NhNumber(3.0, 2.0, 2.0, 0.0), Float64.NaN, NhNumber(Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN)),

    // Handled by math.Pow tests:
    //
    // Pow(±0, y) = ±Inf for y an odd integer < 0
    // Pow(±0, -Inf) = +Inf
    // Pow(±0, +Inf) = +0
    // Pow(±0, y) = +Inf for finite y < 0 and not an odd integer
    // Pow(±0, y) = ±0 for y an odd integer > 0
    // Pow(±0, y) = +0 for finite y > 0 and not an odd integer
    // Pow(-1, ±Inf) = 1

    // PowReal(x+0ϵ₁+0ϵ₂, +Inf) = +Inf+NaNϵ₁+NaNϵ₂+NaNϵ₁ϵ₂ for |x| > 1
    PowRealTest(NhNumber(2.0, 0.0, 0.0, 0.0), Float64.Inf, NhNumber(Float64.Inf, Float64.NaN, Float64.NaN, Float64.NaN)),
    PowRealTest(NhNumber(3.0, 0.0, 0.0, 0.0), Float64.Inf, NhNumber(Float64.Inf, Float64.NaN, Float64.NaN, Float64.NaN)),

    // PowReal(x+xϵ₁+yϵ₂, +Inf) = +Inf+Infϵ₁+Infϵ₂+NaNϵ₁ϵ₂ for |x| > 1
    PowRealTest(NhNumber(2.0, 1.0, 1.0, 0.0), Float64.Inf, NhNumber(Float64.Inf, Float64.Inf, Float64.Inf, Float64.NaN)),
    PowRealTest(NhNumber(3.0, 1.0, 1.0, 0.0), Float64.Inf, NhNumber(Float64.Inf, Float64.Inf, Float64.Inf, Float64.NaN)),
    PowRealTest(NhNumber(2.0, 2.0, 2.0, 0.0), Float64.Inf, NhNumber(Float64.Inf, Float64.Inf, Float64.Inf, Float64.NaN)),
    PowRealTest(NhNumber(3.0, 2.0, 2.0, 0.0), Float64.Inf, NhNumber(Float64.Inf, Float64.Inf, Float64.Inf, Float64.NaN)),
    PowRealTest(NhNumber(3.0, 2.0, 3.0, 0.0), Float64.Inf, NhNumber(Float64.Inf, Float64.Inf, Float64.Inf, Float64.NaN)),

    // PowReal(x, -Inf) = +0+NaNϵ₁+NaNϵ₂+NaNϵ₁ϵ₂ for |x| > 1
    PowRealTest(NhNumber(2.0, 0.0, 0.0, 0.0), -Float64.Inf, NhNumber(0.0, Float64.NaN, Float64.NaN, Float64.NaN)),
    PowRealTest(NhNumber(3.0, 0.0, 0.0, 0.0), -Float64.Inf, NhNumber(0.0, Float64.NaN, Float64.NaN, Float64.NaN)),
    PowRealTest(NhNumber(2.0, 1.0, 1.0, 0.0), -Float64.Inf, NhNumber(0.0, Float64.NaN, Float64.NaN, Float64.NaN)),
    PowRealTest(NhNumber(3.0, 1.0, 1.0, 0.0), -Float64.Inf, NhNumber(0.0, Float64.NaN, Float64.NaN, Float64.NaN)),
    PowRealTest(NhNumber(2.0, 2.0, 2.0, 0.0), -Float64.Inf, NhNumber(0.0, Float64.NaN, Float64.NaN, Float64.NaN)),
    PowRealTest(NhNumber(3.0, 2.0, 2.0, 0.0), -Float64.Inf, NhNumber(0.0, Float64.NaN, Float64.NaN, Float64.NaN)),

    // PowReal(x+yϵ₁+zϵ₂, +Inf) = +0+NaNϵ₁+NaNϵ₂+NaNϵ₁ϵ₂ for |x| < 1
    PowRealTest(NhNumber(0.1, 0.0, 0.0, 0.0), Float64.Inf, NhNumber(0.0, Float64.NaN, Float64.NaN, Float64.NaN)),
    PowRealTest(NhNumber(0.1, 0.1, 0.1, 0.0), Float64.Inf, NhNumber(0.0, Float64.NaN, Float64.NaN, Float64.NaN)),
    PowRealTest(NhNumber(0.2, 0.2, 0.2, 0.0), Float64.Inf, NhNumber(0.0, Float64.NaN, Float64.NaN, Float64.NaN)),
    PowRealTest(NhNumber(0.5, 0.3, 0.5, 0.0), Float64.Inf, NhNumber(0.0, Float64.NaN, Float64.NaN, Float64.NaN)),

    // PowReal(x+0ϵ₁+0ϵ₂, -Inf) = +Inf+NaNϵ₁+NaNϵ₂+NaNϵ₁ϵ₂ for |x| < 1
    PowRealTest(NhNumber(0.1, 0.0, 0.0, 0.0), -Float64.Inf, NhNumber(Float64.Inf, Float64.NaN, Float64.NaN, Float64.NaN)
    ),
    PowRealTest(NhNumber(0.2, 0.0, 0.0, 0.0), -Float64.Inf, NhNumber(Float64.Inf, Float64.NaN, Float64.NaN, Float64.NaN)
    ),

    // PowReal(x, -Inf) = +Inf-Infϵ₁-Infϵ₂+NaNϵ₁ϵ₂ for |x| < 1
    PowRealTest(NhNumber(0.1, 0.1, 0.1, 0.0), -Float64.Inf, NhNumber(Float64.Inf, -Float64.Inf, -Float64.Inf,
        Float64.NaN)),
    PowRealTest(NhNumber(0.2, 0.1, 0.1, 0.0), -Float64.Inf, NhNumber(Float64.Inf, -Float64.Inf, -Float64.Inf,
        Float64.NaN)),
    PowRealTest(NhNumber(0.1, 0.2, 0.2, 0.0), -Float64.Inf, NhNumber(Float64.Inf, -Float64.Inf, -Float64.Inf,
        Float64.NaN)),
    PowRealTest(NhNumber(0.2, 0.3, 0.2, 0.0), -Float64.Inf, NhNumber(Float64.Inf, -Float64.Inf, -Float64.Inf,
        Float64.NaN)),
    PowRealTest(NhNumber(0.1, 1.0, 1.0, 0.0), -Float64.Inf, NhNumber(Float64.Inf, -Float64.Inf, -Float64.Inf,
        Float64.NaN)),
    PowRealTest(NhNumber(0.2, 1.0, 1.0, 0.0), -Float64.Inf, NhNumber(Float64.Inf, -Float64.Inf, -Float64.Inf,
        Float64.NaN)),
    PowRealTest(NhNumber(0.1, 2.0, 2.0, 0.0), -Float64.Inf, NhNumber(Float64.Inf, -Float64.Inf, -Float64.Inf,
        Float64.NaN)),
    PowRealTest(NhNumber(0.2, 2.0, 2.0, 0.0), -Float64.Inf, NhNumber(Float64.Inf, -Float64.Inf, -Float64.Inf,
        Float64.NaN)),

    // Handled by math.Pow tests:
    //
    // Pow(+Inf, y) = +Inf for y > 0
    // Pow(+Inf, y) = +0 for y < 0
    // Pow(-Inf, y) = Pow(-0, -y)

    // PowReal(x, y) = NaN+NaNϵ₁+NaNϵ₂+NaNϵ₁ϵ₂ for finite x < 0 and finite non-integer y
    PowRealTest(NhNumber(-1.0, -1.0, -1.0, 0.0), 0.5, NhNumber(Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN)),
    PowRealTest(NhNumber(-1.0, 2.0, 2.0, 0.0), 0.5, NhNumber(Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN)),
    PowRealTest(NhNumber(-1.0, -1.0, 2.0, 0.0), 0.5, NhNumber(Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN))
]

@Test
public func testPowReal(): Unit {
    let tol = 1e-15
    for (test in powRealTests) {
        let got = nhPowReal(test.d, test.p)
        @Expect(sameHyperdual(got, test.want, tol))
    }
}

private func sameHyperdual(a: NhNumber, b: NhNumber, tol: Float64): Bool {
    return same(a.real, b.real, tol) && same(a.e1mag, b.e1mag, tol) && same(a.e2mag, b.e2mag, tol) && same(a.e1E2mag,
        b.e1E2mag, tol)
}

/*
因为仓颉语言与go语言不同(go语言中-0.0==0.0,仓颉不是)，所以这里额外添加了对0.0,-0.0的判断
 */
private func same(a: Float64, b: Float64, tol: Float64): Bool {
    if ((a == 0.0 && b == -0.0) || (a == -0.0 && b == 0.0)) {
        return true
    }
    return (a.isNaN() && b.isNaN()) || (equalWithinAbsOrRel(a, b, tol, tol)) && (a.toBits() & (1 << 63)) == (b.toBits() &
        (1 << 63))
}
