package cjnum.num.dual

import cjnum.floats.scalar.*
import cjnum.util.*
import std.math.*
import std.unittest.*
import std.unittest.testmacro.*

// TODO var formatTests = []struct

// TODO func TestFormat(t *testing.T) {

/*
因为dAtanh需要使用参数negZero所以negZero前置
 */
private let negZero = Float64.fromBits(1 << 63)

private func dSin(x: Float64): Float64 {
    return cos(x)
}

private func dCos(x: Float64): Float64 {
    return -sin(x)
}

private func dTan(x: Float64): Float64 {
    return sec(x) * sec(x)
}

/*
因为仓颉math库函数sqrt处理小于0的数的方式是直接抛出异常终止程序，
不利于程序运行，所以这里修改了dAsin的逻辑能对输入值进行判断，
如果sqrt的输入值小于0就返回NaN，否则正常使用math库函数sqrt正常处理范围内的数
 */
private func dAsin(x: Float64): Float64 {
    if (1.0 - x * x < 0.0) {
        return Float64.NaN
    }
    return 1.0 / sqrt(1.0 - x * x)
}

/*
因为仓颉math库函数sqrt处理小于0的数的方式是直接抛出异常终止程序，
不利于程序运行，所以这里修改了dAcos的逻辑能对输入值进行判断，
如果sqrt的输入值小于0就返回NaN，否则正常使用math库函数sqrt正常处理范围内的数
 */
private func dAcos(x: Float64): Float64 {
    if (1.0 - x * x < 0.0) {
        return Float64.NaN
    }
    return -1.0 / sqrt(1.0 - x * x)
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

/*
因为仓颉math库函数sqrt处理小于0的数的方式是直接抛出异常终止程序，
不利于程序运行，所以这里修改了dAcosh的逻辑能对输入值进行判断，
如果sqrt的输入值x-1小于0就返回NaN，否则正常使用math库函数sqrt正常处理范围内的数
 */
private func dAcosh(x: Float64): Float64 {
    if (x - 1.0 < 0.0) {
        return Float64.NaN
    }
    return 1.0 / (sqrt(x - 1.0) * sqrt(x + 1.0))
}

private func dAtanh(x: Float64): Float64 {
    if (abs(x) == 1.0) {
        return Float64.NaN
    } else if (x.isInf()) { //同时检查正负无穷大
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

/*
因为仓颉math库函数sqrt处理小于0的数的方式是直接抛出异常终止程序，
不利于程序运行，所以这里修改了dASqrt的逻辑把sqrt换成了nanSqrt，
如果nanSqrt的输入值x小于0就返回NaN，否则正常使用math库函数sqrt正常处理范围内的数
 */
private func dSqrt(x: Float64): Float64 {
    // For whatever reason, math.Sqrt(-0) returns -0.
    // In this case, that is clearly a wrong approach.
    if (x == 0.0) {
        return Float64.Inf
    }
    return 0.5 / nanSqrt(x)
}

private func dInv(x: Float64): Float64 {
    return -1.0 / (x * x)
}

// Helpers:

private func sec(x: Float64): Float64 {
    return 1.0 / cos(x)
}

private func sech(x: Float64): Float64 {
    return 1.0 / cosh(x)
}

private class DualTest {
    let name: String
    let x: Array<Float64>
    let fnDual: (x: NdNumber) -> NdNumber
    let fn: (x: Float64) -> Float64
    let dFn: (x: Float64) -> Float64
    public init(name: String, x: Array<Float64>, fnDual: (x: NdNumber) -> NdNumber, fn: (x: Float64) -> Float64,
        dFn: (x: Float64) -> Float64) {
        this.name = name
        this.x = x
        this.fnDual = fnDual
        this.fn = fn
        this.dFn = dFn
    }
}

private let dualTests: Array<DualTest> = [
    DualTest(
        "sin",
        [Float64.NaN, -Float64.Inf, -3.0, -2.0, -1.0, -0.5, negZero, 0.0, 0.5, 1.0, 2.0, 3.0, Float64.Inf],
        ndSin,
        sin,
        dSin
    ),
    DualTest(
        "cos",
        [Float64.NaN, -Float64.Inf, -3.0, -2.0, -1.0, -0.5, negZero, 0.0, 0.5, 1.0, 2.0, 3.0, Float64.Inf],
        ndCos,
        cos,
        dCos
    ),
    DualTest(
        "tan",
        [Float64.NaN, -Float64.Inf, -3.0, -2.0, -1.0, -0.5, negZero, 0.0, 0.5, 1.0, 2.0, 3.0, Float64.Inf],
        ndTan,
        tan,
        dTan
    ),
    DualTest(
        "sinh",
        [Float64.NaN, -Float64.Inf, -3.0, -2.0, -1.0, -0.5, negZero, 0.0, 0.5, 1.0, 2.0, 3.0, Float64.Inf],
        ndSinh,
        sinh,
        dSinh
    ),
    DualTest(
        "cosh",
        [Float64.NaN, -Float64.Inf, -3.0, -2.0, -1.0, -0.5, negZero, 0.0, 0.5, 1.0, 2.0, 3.0, Float64.Inf],
        ndCosh,
        cosh,
        dCosh
    ),
    DualTest(
        "tanh",
        [Float64.NaN, -Float64.Inf, -3.0, -2.0, -1.0, -0.5, negZero, 0.0, 0.5, 1.0, 2.0, 3.0, Float64.Inf],
        ndTanh,
        tanh,
        dTanh
    ),
    DualTest(
        "asin",
        [Float64.NaN, -Float64.Inf, -3.0, -2.0, -1.0, -0.5, negZero, 0.0, 0.5, 1.0, 2.0, 3.0, Float64.Inf],
        ndAsin,
        nanAsin,
        dAsin
    ),
    DualTest(
        "acos",
        [Float64.NaN, -Float64.Inf, -3.0, -2.0, -1.0, -0.5, negZero, 0.0, 0.5, 1.0, 2.0, 3.0, Float64.Inf],
        ndAcos,
        nanAcos,
        dAcos
    ),
    DualTest(
        "atan",
        [Float64.NaN, -Float64.Inf, -3.0, -2.0, -1.0, -0.5, negZero, 0.0, 0.5, 1.0, 2.0, 3.0, Float64.Inf],
        ndAtan,
        atan,
        dAtan
    ),
    DualTest(
        "asinh",
        [Float64.NaN, -Float64.Inf, -3.0, -2.0, -1.0, -0.5, negZero, 0.0, 0.5, 1.0, 2.0, 3.0, Float64.Inf],
        ndAsinh,
        asinh,
        dAsinh
    ),
    DualTest(
        "acosh",
        [Float64.NaN, -Float64.Inf, -3.0, -2.0, -1.0, -0.5, negZero, 0.0, 0.5, 1.0, 2.0, 3.0, Float64.Inf],
        ndAcosh,
        nanAcosh,
        dAcosh
    ),
    DualTest(
        "atanh",
        [Float64.NaN, -Float64.Inf, -3.0, -2.0, -1.0, -0.5, negZero, 0.0, 0.5, 1.0, 2.0, 3.0, Float64.Inf],
        ndAtanh,
        nanAtanh,
        dAtanh
    ),
    DualTest(
        "exp",
        [Float64.NaN, -Float64.Inf, -3.0, -2.0, -1.0, -0.5, negZero, 0.0, 0.5, 1.0, 2.0, 3.0, Float64.Inf],
        ndExp,
        exp,
        dExp
    ),
    DualTest(
        "log",
        [Float64.NaN, -Float64.Inf, -3.0, -2.0, -1.0, -0.5, negZero, 0.0, 0.5, 1.0, 2.0, 3.0, Float64.Inf],
        ndLog,
        log,
        dLog
    ),
    DualTest(
        "inv",
        [Float64.NaN, -Float64.Inf, -3.0, -2.0, -1.0, -0.5, negZero, 0.0, 0.5, 1.0, 2.0, 3.0, Float64.Inf],
        ndInv,
        {x => 1.0 / x},
        dInv,
    ),
    DualTest(
        "sqrt",
        [Float64.NaN, -Float64.Inf, -3.0, -2.0, -1.0, -0.5, negZero, 0.0, 0.5, 1.0, 2.0, 3.0, Float64.Inf],
        ndSqrt,
        nanSqrt,
        dSqrt,
    ),
    DualTest(
        "Fike example fn",
        [1.0, 2.0, 3.0, 4.0, 5.0],
        {
            x => ndMul(
                ndExp(x),
                ndInv(
                    ndSqrt(
                        ndAdd(
                            ndPowReal(ndSin(x), 3.0),
                            ndPowReal(ndCos(x), 3.0)
                        )
                    )
                )
            )
        },
        {x => exp(x) / nanSqrt(pow(sin(x), 3.0) + pow(cos(x), 3.0))},
        {
            x => exp(x) * (3.0 * cos(x) + 5.0 * cos(3.0 * x) + 9.0 * sin(x) + sin(3.0 * x)) / (8.0 * pow(
                pow(sin(x), 3.0) + pow(cos(x), 3.0), 1.5))
        },
    )
]

@Test
public func testDual(): Unit {
    let tol = 1e-15
    for (dualTest in dualTests) {
        for (x in dualTest.x) {
            let fxDual = dualTest.fnDual(NdNumber(x, 1.0))
            let fx = dualTest.fn(x)
            let dFx = dualTest.dFn(x)
            @Expect(same(fxDual.real, fx, tol))
            @Expect(same(fxDual.emag, dFx, tol))
        }
    }
}

private class PowRealTest {
    let d: NdNumber
    let p: Float64
    let want: NdNumber
    public init(d: NdNumber, p: Float64, want: NdNumber) {
        this.d = d
        this.p = p
        this.want = want
    }
}

private let powRealTests: Array<PowRealTest> = [
    // PowReal(NaN+xϵ, ±0) = 1+NaNϵ for any x
    PowRealTest(
        NdNumber(Float64.NaN, 0.0),
        0.0,
        NdNumber(1.0, Float64.NaN)
    ),
    PowRealTest(
        NdNumber(Float64.NaN, 0.0),
        negZero,
        NdNumber(1.0, Float64.NaN)
    ),
    PowRealTest(
        NdNumber(Float64.NaN, 1.0),
        0.0,
        NdNumber(1.0, Float64.NaN)
    ),
    PowRealTest(
        NdNumber(Float64.NaN, 2.0),
        negZero,
        NdNumber(1.0, Float64.NaN)
    ),
    PowRealTest(
        NdNumber(Float64.NaN, 3.0),
        0.0,
        NdNumber(1.0, Float64.NaN)
    ),
    PowRealTest(
        NdNumber(Float64.NaN, 1.0),
        negZero,
        NdNumber(1.0, Float64.NaN)
    ),
    PowRealTest(
        NdNumber(Float64.NaN, 2.0),
        0.0,
        NdNumber(1.0, Float64.NaN)
    ),
    PowRealTest(
        NdNumber(Float64.NaN, 3.0),
        negZero,
        NdNumber(1.0, Float64.NaN)
    ),

    // PowReal(x, ±0) = 1 for any x
    PowRealTest(
        NdNumber(0.0, 0.0),
        0.0,
        NdNumber(1.0, 0.0)
    ),
    PowRealTest(
        NdNumber(negZero, 0.0),
        negZero,
        NdNumber(1.0, 0.0)
    ),
    PowRealTest(
        NdNumber(Float64.Inf, 0.0),
        0.0,
        NdNumber(1.0, 0.0)
    ),
    PowRealTest(
        NdNumber(-Float64.Inf, 0.0),
        negZero,
        NdNumber(1.0, 0.0)
    ),
    PowRealTest(
        NdNumber(0.0, 1.0),
        0.0,
        NdNumber(1.0, 0.0)
    ),
    PowRealTest(
        NdNumber(negZero, 1.0),
        negZero,
        NdNumber(1.0, 0.0)
    ),
    PowRealTest(
        NdNumber(Float64.Inf, 1.0),
        0.0,
        NdNumber(1.0, 0.0)
    ),
    PowRealTest(
        NdNumber(-Float64.Inf, 1.0),
        negZero,
        NdNumber(1.0, 0.0)
    ),

    // PowReal(1+xϵ, y) = (1+xyϵ) for any y
    PowRealTest(
        NdNumber(1.0, 0.0),
        0.0,
        NdNumber(1.0, 0.0)
    ),
    PowRealTest(
        NdNumber(1.0, 0.0),
        1.0,
        NdNumber(1.0, 0.0)
    ),
    PowRealTest(
        NdNumber(1.0, 0.0),
        2.0,
        NdNumber(1.0, 0.0)
    ),
    PowRealTest(
        NdNumber(1.0, 0.0),
        3.0,
        NdNumber(1.0, 0.0)
    ),
    PowRealTest(
        NdNumber(1.0, 1.0),
        0.0,
        NdNumber(1.0, 0.0)
    ),
    PowRealTest(
        NdNumber(1.0, 1.0),
        1.0,
        NdNumber(1.0, 1.0)
    ),
    PowRealTest(
        NdNumber(1.0, 1.0),
        2.0,
        NdNumber(1.0, 2.0)
    ),
    PowRealTest(
        NdNumber(1.0, 1.0),
        3.0,
        NdNumber(1.0, 3.0)
    ),
    PowRealTest(
        NdNumber(1.0, 2.0),
        0.0,
        NdNumber(1.0, 0.0)
    ),
    PowRealTest(
        NdNumber(1.0, 2.0),
        1.0,
        NdNumber(1.0, 2.0)
    ),
    PowRealTest(
        NdNumber(1.0, 2.0),
        2.0,
        NdNumber(1.0, 4.0)
    ),
    PowRealTest(
        NdNumber(1.0, 2.0),
        3.0,
        NdNumber(1.0, 6.0)
    ),

    // PowReal(x, 1) = x for any x
    PowRealTest(
        NdNumber(0.0, 0.0),
        1.0,
        NdNumber(0.0, 0.0)
    ),
    PowRealTest(
        NdNumber(negZero, 0.0),
        1.0,
        NdNumber(negZero, 0.0)
    ),
    PowRealTest(
        NdNumber(0.0, 1.0),
        1.0,
        NdNumber(0.0, 1.0)
    ),
    PowRealTest(
        NdNumber(negZero, 1.0),
        1.0,
        NdNumber(negZero, 1.0)
    ),
    PowRealTest(
        NdNumber(Float64.NaN, 0.0),
        1.0,
        NdNumber(Float64.NaN, 0.0)
    ),
    PowRealTest(
        NdNumber(Float64.NaN, 1.0),
        1.0,
        NdNumber(Float64.NaN, 1.0)
    ),
    PowRealTest(
        NdNumber(Float64.NaN, 2.0),
        1.0,
        NdNumber(Float64.NaN, 2.0)
    ),

    // PowReal(NaN+xϵ, y) = NaN+NaNϵ
    PowRealTest(
        NdNumber(Float64.NaN, 0.0),
        2.0,
        NdNumber(Float64.NaN, Float64.NaN)
    ),
    PowRealTest(
        NdNumber(Float64.NaN, 0.0),
        3.0,
        NdNumber(Float64.NaN, Float64.NaN)
    ),
    PowRealTest(
        NdNumber(Float64.NaN, 1.0),
        2.0,
        NdNumber(Float64.NaN, Float64.NaN)
    ),
    PowRealTest(
        NdNumber(Float64.NaN, 1.0),
        3.0,
        NdNumber(Float64.NaN, Float64.NaN)
    ),
    PowRealTest(
        NdNumber(Float64.NaN, 2.0),
        2.0,
        NdNumber(Float64.NaN, Float64.NaN)
    ),
    PowRealTest(
        NdNumber(Float64.NaN, 2.0),
        3.0,
        NdNumber(Float64.NaN, Float64.NaN)
    ),

    // PowReal(x, NaN) = NaN+NaNϵ
    PowRealTest(
        NdNumber(0.0, 0.0),
        Float64.NaN,
        NdNumber(Float64.NaN, Float64.NaN)
    ),
    PowRealTest(
        NdNumber(2.0, 0.0),
        Float64.NaN,
        NdNumber(Float64.NaN, Float64.NaN)
    ),
    PowRealTest(
        NdNumber(3.0, 0.0),
        Float64.NaN,
        NdNumber(Float64.NaN, Float64.NaN)
    ),
    PowRealTest(
        NdNumber(0.0, 1.0),
        Float64.NaN,
        NdNumber(Float64.NaN, Float64.NaN)
    ),
    PowRealTest(
        NdNumber(2.0, 1.0),
        Float64.NaN,
        NdNumber(Float64.NaN, Float64.NaN)
    ),
    PowRealTest(
        NdNumber(3.0, 1.0),
        Float64.NaN,
        NdNumber(Float64.NaN, Float64.NaN)
    ),
    PowRealTest(
        NdNumber(0.0, 2.0),
        Float64.NaN,
        NdNumber(Float64.NaN, Float64.NaN)
    ),
    PowRealTest(
        NdNumber(2.0, 2.0),
        Float64.NaN,
        NdNumber(Float64.NaN, Float64.NaN)
    ),
    PowRealTest(
        NdNumber(3.0, 2.0),
        Float64.NaN,
        NdNumber(Float64.NaN, Float64.NaN)
    ),

    // Handled by math.Pow tests:
    //
    // Pow(±0, y) = ±Inf for y an odd integer < 0
    // Pow(±0, -Inf) = +Inf
    // Pow(±0, +Inf) = +0
    // Pow(±0, y) = +Inf for finite y < 0 and not an odd integer
    // Pow(±0, y) = ±0 for y an odd integer > 0
    // Pow(±0, y) = +0 for finite y > 0 and not an odd integer
    // Pow(-1, ±Inf) = 1 

    // PowReal(x+0ϵ, +Inf) = +Inf+NaNϵ for |x| > 1
    PowRealTest(NdNumber(2.0, 0.0), Float64.Inf, NdNumber(Float64.Inf, Float64.NaN)),
    PowRealTest(NdNumber(3.0, 0.0), Float64.Inf, NdNumber(Float64.Inf, Float64.NaN)),

    // PowReal(x+yϵ, +Inf) = +Inf for |x| > 1
    PowRealTest(NdNumber(2.0, 1.0), Float64.Inf, NdNumber(Float64.Inf, Float64.Inf)),
    PowRealTest(NdNumber(3.0, 1.0), Float64.Inf, NdNumber(Float64.Inf, Float64.Inf)),
    PowRealTest(NdNumber(2.0, 2.0), Float64.Inf, NdNumber(Float64.Inf, Float64.Inf)),
    PowRealTest(NdNumber(3.0, 2.0), Float64.Inf, NdNumber(Float64.Inf, Float64.Inf)),

    // PowReal(x, -Inf) = +0+NaNϵ for |x| > 1
    PowRealTest(NdNumber(2.0, 0.0), -Float64.Inf, NdNumber(0.0, Float64.NaN)),
    PowRealTest(NdNumber(3.0, 0.0), -Float64.Inf, NdNumber(0.0, Float64.NaN)),
    PowRealTest(NdNumber(2.0, 1.0), -Float64.Inf, NdNumber(0.0, Float64.NaN)),
    PowRealTest(NdNumber(3.0, 1.0), -Float64.Inf, NdNumber(0.0, Float64.NaN)),
    PowRealTest(NdNumber(2.0, 2.0), -Float64.Inf, NdNumber(0.0, Float64.NaN)),
    PowRealTest(NdNumber(3.0, 2.0), -Float64.Inf, NdNumber(0.0, Float64.NaN)),

    // PowReal(x+yϵ, +Inf) = +0+NaNϵ for |x| < 1
    PowRealTest(NdNumber(0.1, 0.0), Float64.Inf, NdNumber(0.0, Float64.NaN)),
    PowRealTest(NdNumber(0.1, 0.1), Float64.Inf, NdNumber(0.0, Float64.NaN)),
    PowRealTest(NdNumber(0.2, 0.2), Float64.Inf, NdNumber(0.0, Float64.NaN)),
    PowRealTest(NdNumber(0.5, 0.5), Float64.Inf, NdNumber(0.0, Float64.NaN)),

    // PowReal(x+0ϵ, -Inf) = +Inf+NaNϵ for |x| < 1
    PowRealTest(NdNumber(0.1, 0.0), -Float64.Inf, NdNumber(Float64.Inf, Float64.NaN)),
    PowRealTest(NdNumber(0.2, 0.0), -Float64.Inf, NdNumber(Float64.Inf, Float64.NaN)),

    // PowReal(x, -Inf) = +Inf-Infϵ for |x| < 1
    PowRealTest(NdNumber(0.1, 0.1), -Float64.Inf, NdNumber(Float64.Inf, -Float64.Inf)),
    PowRealTest(NdNumber(0.2, 0.1), -Float64.Inf, NdNumber(Float64.Inf, -Float64.Inf)),
    PowRealTest(NdNumber(0.1, 0.2), -Float64.Inf, NdNumber(Float64.Inf, -Float64.Inf)),
    PowRealTest(NdNumber(0.2, 0.2), -Float64.Inf, NdNumber(Float64.Inf, -Float64.Inf)),
    PowRealTest(NdNumber(0.1, 1.0), -Float64.Inf, NdNumber(Float64.Inf, -Float64.Inf)),
    PowRealTest(NdNumber(0.2, 1.0), -Float64.Inf, NdNumber(Float64.Inf, -Float64.Inf)),
    PowRealTest(NdNumber(0.1, 2.0), -Float64.Inf, NdNumber(Float64.Inf, -Float64.Inf)),
    PowRealTest(NdNumber(0.2, 2.0), -Float64.Inf, NdNumber(Float64.Inf, -Float64.Inf)),

    // Handled by math.Pow tests:
    //
    // Pow(+Inf, y) = +Inf for y > 0
    // Pow(+Inf, y) = +0 for y < 0
    // Pow(-Inf, y) = Pow(-0, -y)

    // PowReal(x, y) = NaN+NaNϵ for finite x < 0 and finite non-integer y
    PowRealTest(NdNumber(-1.0, -1.0), 0.5, NdNumber(Float64.NaN, Float64.NaN)),
    PowRealTest(NdNumber(-1.0, 2.0), 0.5, NdNumber(Float64.NaN, Float64.NaN))
]

@Test
public func testPowReal(): Unit {
    let tol = 1e-15
    for (test in powRealTests) {
        let got = ndPowReal(test.d, test.p)
        @Expect(sameDual(got, test.want, tol))
    }
}

private func sameDual(a: NdNumber, b: NdNumber, tol: Float64): Bool {
    return same(a.real, b.real, tol) && same(a.emag, b.emag, tol)
}

private func same(a: Float64, b: Float64, tol: Float64): Bool {
    return (a.isNaN() && b.isNaN()) || equalWithinAbsOrRel(a, b, tol, tol)
}

// @Test
// @Parallel
// public func testScale(): Unit {
//     for (test in [
//         (2.0, NdNumber(1.0, 3.0), NdNumber(2.0, 6.0)),
//         (-1.0, NdNumber(1.0, -2.0), NdNumber(-1.0, 2.0)),
//         (0.0, NdNumber(5.0, 7.0), NdNumber(0.0, 0.0)),
//         (1.0, NdNumber(-3.0, 4.0), NdNumber(-3.0, 4.0))
//     ]) {
//         let got = ndScale(test[0], test[1])
//         if (!sameDual(got, test[2], 1e-15)) {
//             throw Exception("Scale test failed")
//         }
//     }
// }

// @Test
// @Parallel
// public func testAbs() {
//     for (test in [
//         (NdNumber(3.0, 2.0), NdNumber(3.0, 2.0)),
//         (NdNumber(-3.0, 2.0), NdNumber(3.0, -2.0)),
//         (NdNumber(0.0, 5.0), NdNumber(0.0, 5.0)),
//         (NdNumber(-0.0, 7.0), NdNumber(-0.0, 7.0))
//     ]) {
//         let got = ndAbs(test[0])
//         if (!(same(got.real, test[1].real, 1e-15) && same(got.emag, test[1].emag, 1e-15))) {
//             throw Exception("Abs test failed")
//         }
//     }
// }