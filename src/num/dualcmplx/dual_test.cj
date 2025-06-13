package cjnum.num.dualcmplx

import cjnum.floats.scalar.*
import cjnum.complex.*
import cjnum.util.*
import std.math.*
import std.unittest.*
import std.unittest.testmacro.*

// TODO var formatTests = []struct {

// TODO func TestFormat(t *testing.T) {

/*
以下的TODO部分因为缺乏相关的测试代码支持所以没有实现
*/
// TODO func sqrt(x complex128) complex128 {

// First derivatives:

private func dExp(x: Complex128): Complex128 {
    if (x.imag == 0.0) {
        return c128Exp(x)
    }
    return (c128Exp(x) - c128Exp(x.conj())) / (x - x.conj())
}

private func dLog(x: Complex128): Complex128 {
    if (x.isInf()) {
        return Complex128(0.0, 0.0)
    }
    if (x.real == 0.0 && x.imag == 0.0) {
        if (copysign(1.0, x.real) < 0.0) {
            return Complex128(-Float64.Inf, Float64.NaN)
        }
        return Complex128(Float64.Inf, Float64.NaN)
    }
    return (c128Log(x) - c128Log(x.conj())) / (x - x.conj())
}

private func dInv(x: Complex128): Complex128 {
    return Complex128(-1.0, 0.0) / (x * x.conj())
}

private let negZero = copysign(0.0, -1.0)
private let zeroCmplx = Complex128(0.0, 0.0)
private let negZeroCmplx = zeroCmplx * (-1.0)
private let one = Complex128(1.0, 1.0)
private let negOne = Complex128(-1.0, -1.0)
private let half = one / 2.0
private let negHalf = negOne / 2.0
private let two = Complex128(2.0, 2.0)
private let negTwo = Complex128(-2.0, -2.0)
private let three = Complex128(3.0, 3.0)
private let negThree = Complex128(-3.0, 3.0)

private class DualTest {
    let name: String
    let x: Array<Complex128>
    let fnDual: (x: NdcNumber) -> NdcNumber
    let fn: (x: Complex128) -> Complex128
    let dFn: (x: Complex128) -> Complex128
    public init(name: String, x: Array<Complex128>, fnDual: (x: NdcNumber) -> NdcNumber, fn: (x: Complex128) -> Complex128,
        dFn: (x: Complex128) -> Complex128) {
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
        [Complex128.NaN, Complex128.Inf, negThree, negTwo, negOne, negHalf, negZeroCmplx, zeroCmplx, half, one, two,
            three],
        ndcExp,
        c128Exp,
        dExp,
    ),
    DualTest(
        "log",
        [Complex128.NaN, Complex128.Inf, negThree, negTwo, negOne, negHalf, negZeroCmplx, zeroCmplx, half, one, two,
            three],
        ndcLog,
        c128Log,
        dLog,
    ),
    DualTest(
        "inv",
        [Complex128.NaN, Complex128.Inf, negThree, negTwo, negOne, negHalf, negZeroCmplx, zeroCmplx, half, one, two,
            three],
        ndcInv,
        {x: Complex128 => Complex128(1.0, 0.0) / x},
        dInv,
    )
    /*
    因为仓颉语言不能存在空的元素，而源项目这段测试用例有函数没有实现（TODO(kortschak): Find a concise dSqrt.），所以这里删除了这段测试用例
     */
    // ,
    // DualTest(
    // 	"sqrt",
    // 	[Complex128.NaN, Complex128.Inf, negThree, negTwo, negOne, negHalf, negZeroCmplx, zeroCmplx, half, one, two, three],
    // 	Sqrt,
    // 	sqrt,
    // 	// TODO(kortschak): Find a concise dSqrt.
    // )
]

@Test
private func testNumber(): Unit {
    let tol = 1e-15
    for (test in dualTests) {
        for (x in test.x) {
            let fxDual = test.fnDual(NdcNumber(x, Complex128(1.0, 0.0)))
            let fx = test.fn(x)
            @Expect(same(fxDual.real, fx, tol))
            let dFx = test.dFn(x)
            @Expect(same(fxDual.dual, dFx, tol))
        }
    }
}

private var invTests: Array<NdcNumber> = [
    NdcNumber(
        Complex128(1.0, 0.0),
        Complex128(0.0, 0.0)
    ),
    NdcNumber(
        Complex128(1.0, 0.0),
        Complex128(1.0, 0.0)
    ),
    NdcNumber(
        Complex128(0.0, 1.0),
        Complex128(1.0, 0.0)
    ),
    NdcNumber(
        Complex128(1.0, 0.0),
        Complex128(1.0, 1.0)
    ),
    NdcNumber(
        Complex128(1.0, 1.0),
        Complex128(1.0, 1.0)
    ),
    NdcNumber(
        Complex128(1.0, 10.0),
        Complex128(1.0, 5.0)
    ),
    NdcNumber(
        Complex128(10.0, 1.0),
        Complex128(5.0, 1.0)
    )
]

@Test
private func testInv(): Unit {
    let tol = 1e-15
    for (x in invTests) {
        let got = ndcMul(x, ndcInv(x))
        let want = NdcNumber(Complex128(1.0, 0.0), Complex128(0.0, 0.0))
        @Expect(sameDual(got, want, tol))
    }
}

private var expLogTests: Array<NdcNumber> = [
    NdcNumber(Complex128(0.0, 1.0), Complex128(0.0, 1.0)),
    NdcNumber(Complex128(1.0, 1.0), Complex128(1.0, 1.0)),
    NdcNumber(Complex128(1.0, 0.1), Complex128(1.0, 1.0)),
    NdcNumber(Complex128(1.0, 0.01), Complex128(1.0, 1.0)),
    NdcNumber(Complex128(1.0, 0.0001), Complex128(1.0, 1.0)),
    NdcNumber(Complex128(1.0, 0.000001), Complex128(1.0, 1.0)),
    NdcNumber(Complex128(1.0, 0.00000001), Complex128(1.0, 1.0)),
    NdcNumber(Complex128(1.0, 0.0000000001), Complex128(1.0, 1.0)),
    NdcNumber(Complex128(1.0, 0.000000000001), Complex128(1.0, 1.0)),
    NdcNumber(Complex128(1.0, 0.00000000000001), Complex128(1.0, 1.0)),
    NdcNumber(Complex128(0.0, 0.0), Complex128(1.0, 1.0)),
    NdcNumber(Complex128(0.0, 0.0), Complex128(1.0, 2.0)),
    NdcNumber(Complex128(0.0, 0.0), Complex128(2.0, 1.0)),
    NdcNumber(Complex128(0.0, 0.0), Complex128(2.0, 2.0)),
    NdcNumber(Complex128(0.0, 0.0), Complex128(1.0, 1.0)),
    NdcNumber(Complex128(0.0, 0.0), Complex128(1.0, 5.0)),
    NdcNumber(Complex128(0.0, 0.0), Complex128(5.0, 1.0)),
    NdcNumber(Complex128(1.0, 0.0), Complex128(1.0, 1.0)),
    NdcNumber(Complex128(1.0, 0.0), Complex128(1.0, 2.0)),
    NdcNumber(Complex128(1.0, 0.0), Complex128(2.0, 1.0)),
    NdcNumber(Complex128(1.0, 0.0), Complex128(2.0, 2.0)),
    NdcNumber(Complex128(1.0, 1.0), Complex128(1.0, 1.0)),
    NdcNumber(Complex128(1.0, 3.0), Complex128(1.0, 5.0)),
    NdcNumber(Complex128(2.0, 1.0), Complex128(5.0, 1.0))
]

@Test
private func testExpLog(): Unit {
    let tol = 1e-15
    for (x in expLogTests) {
        let got = ndcLog(ndcExp(x))
        let want = x
        @Expect(sameDual(got, want, tol))
    }
}

@Test
private func testLogExp(): Unit {
    let tol = 1e-15
    for (x in expLogTests) {
        if (x.real == Complex128(0.0, 0.0)) {
            continue
        }
        let got = ndcExp(ndcLog(x))
        let want = x
        @Expect(sameDual(got, want, tol))
    }
}

private class PowRealSpecialTest {
    let d: NdcNumber
    let p: Float64
    let want: NdcNumber
    public init(d: NdcNumber, p: Float64, want: NdcNumber) {
        this.d = d
        this.p = p
        this.want = want
    }
}

private let powRealSpecialTests: Array<PowRealSpecialTest> = [
    // PowReal(NaN+xϵ, ±0) = 1+NaNϵ for any x
    PowRealSpecialTest(NdcNumber(Complex128.NaN, Complex128(0.0, 0.0)), 0.0, NdcNumber(Complex128(1.0, 0.0), Complex128.NaN)),
    PowRealSpecialTest(NdcNumber(Complex128.NaN, Complex128(0.0, 0.0)), negZero, NdcNumber(Complex128(1.0, 0.0),
        Complex128.NaN)),
    PowRealSpecialTest(NdcNumber(Complex128.NaN, Complex128(1.0, 0.0)), 0.0, NdcNumber(Complex128(1.0, 0.0), Complex128.NaN)),
    PowRealSpecialTest(NdcNumber(Complex128.NaN, Complex128(2.0, 0.0)), negZero, NdcNumber(Complex128(1.0, 0.0),
        Complex128.NaN)),
    PowRealSpecialTest(NdcNumber(Complex128.NaN, Complex128(3.0, 0.0)), 0.0, NdcNumber(Complex128(1.0, 0.0), Complex128.NaN)),
    PowRealSpecialTest(NdcNumber(Complex128.NaN, Complex128(1.0, 0.0)), negZero, NdcNumber(Complex128(1.0, 0.0),
        Complex128.NaN)),
    PowRealSpecialTest(NdcNumber(Complex128.NaN, Complex128(2.0, 0.0)), 0.0, NdcNumber(Complex128(1.0, 0.0), Complex128.NaN)),
    PowRealSpecialTest(NdcNumber(Complex128.NaN, Complex128(3.0, 0.0)), negZero, NdcNumber(Complex128(1.0, 0.0),
        Complex128.NaN)),

    // Pow(0+xϵ, y) = 0+Infϵ for all y < 1.
    PowRealSpecialTest(NdcNumber(Complex128(0.0, 0.0), Complex128(0.0, 0.0)), 0.1,
        NdcNumber(Complex128(0.0, 0.0), Complex128.Inf)),
    PowRealSpecialTest(NdcNumber(Complex128(0.0, 0.0), Complex128(0.0, 0.0)), -1.0,
        NdcNumber(Complex128(0.0, 0.0), Complex128.Inf)),
    PowRealSpecialTest(NdcNumber(Complex128(0.0, 0.0), Complex128(1.0, 0.0)), 0.1,
        NdcNumber(Complex128(0.0, 0.0), Complex128.Inf)),
    PowRealSpecialTest(NdcNumber(Complex128(0.0, 0.0), Complex128(1.0, 0.0)), -1.0,
        NdcNumber(Complex128(0.0, 0.0), Complex128.Inf)),
    PowRealSpecialTest(NdcNumber(Complex128(0.0, 0.0), Complex128(1.0, 1.0)), 0.1,
        NdcNumber(Complex128(0.0, 0.0), Complex128.Inf)),
    PowRealSpecialTest(NdcNumber(Complex128(0.0, 0.0), Complex128(1.0, 1.0)), -1.0,
        NdcNumber(Complex128(0.0, 0.0), Complex128.Inf)),
    PowRealSpecialTest(NdcNumber(Complex128(0.0, 0.0), Complex128(0.0, 1.0)), 0.1,
        NdcNumber(Complex128(0.0, 0.0), Complex128.Inf)),
    PowRealSpecialTest(NdcNumber(Complex128(0.0, 0.0), Complex128(0.0, 1.0)), -1.0,
        NdcNumber(Complex128(0.0, 0.0), Complex128.Inf)),
    // Pow(0+xϵ, y) = 0 for all y > 1.
    PowRealSpecialTest(NdcNumber(Complex128(0.0, 0.0), Complex128(0.0, 0.0)), 1.1,
        NdcNumber(Complex128(0.0, 0.0), Complex128(0.0, 0.0))),
    PowRealSpecialTest(NdcNumber(Complex128(0.0, 0.0), Complex128(0.0, 0.0)), 2.0,
        NdcNumber(Complex128(0.0, 0.0), Complex128(0.0, 0.0))),
    PowRealSpecialTest(NdcNumber(Complex128(0.0, 0.0), Complex128(1.0, 0.0)), 1.1,
        NdcNumber(Complex128(0.0, 0.0), Complex128(0.0, 0.0))),
    PowRealSpecialTest(NdcNumber(Complex128(0.0, 0.0), Complex128(1.0, 0.0)), 2.0,
        NdcNumber(Complex128(0.0, 0.0), Complex128(0.0, 0.0))),
    PowRealSpecialTest(NdcNumber(Complex128(0.0, 0.0), Complex128(1.0, 1.0)), 1.1,
        NdcNumber(Complex128(0.0, 0.0), Complex128(0.0, 0.0))),
    PowRealSpecialTest(NdcNumber(Complex128(0.0, 0.0), Complex128(1.0, 1.0)), 2.0,
        NdcNumber(Complex128(0.0, 0.0), Complex128(0.0, 0.0))),
    PowRealSpecialTest(NdcNumber(Complex128(0.0, 0.0), Complex128(0.0, 1.0)), 1.1,
        NdcNumber(Complex128(0.0, 0.0), Complex128(0.0, 0.0))),
    PowRealSpecialTest(NdcNumber(Complex128(0.0, 0.0), Complex128(0.0, 1.0)), 2.0,
        NdcNumber(Complex128(0.0, 0.0), Complex128(0.0, 0.0))),

    // PowReal(x, ±0) = 1 for any x
    PowRealSpecialTest(NdcNumber(Complex128(0.0, 0.0), Complex128(0.0, 0.0)), 0.0,
        NdcNumber(Complex128(1.0, 0.0), Complex128(0.0, 0.0))),
    PowRealSpecialTest(NdcNumber(negZeroCmplx, Complex128(0.0, 0.0)), negZero,
        NdcNumber(Complex128(1.0, 0.0), Complex128(0.0, 0.0))),
    PowRealSpecialTest(NdcNumber(Complex128(Float64.Inf, 0.0), Complex128(0.0, 0.0)), 0.0,
        NdcNumber(Complex128(1.0, 0.0), Complex128(0.0, 0.0))),
    PowRealSpecialTest(NdcNumber(Complex128(Float64.Inf, 0.0), Complex128(0.0, 0.0)), negZero,
        NdcNumber(Complex128(1.0, 0.0), Complex128(0.0, 0.0))),
    PowRealSpecialTest(NdcNumber(Complex128(0.0, 0.0), Complex128(1.0, 0.0)), 0.0,
        NdcNumber(Complex128(1.0, 0.0), Complex128(0.0, 0.0))),
    PowRealSpecialTest(NdcNumber(negZeroCmplx, Complex128(1.0, 0.0)), negZero,
        NdcNumber(Complex128(1.0, 0.0), Complex128(0.0, 0.0))),
    PowRealSpecialTest(NdcNumber(Complex128(Float64.Inf, 0.0), Complex128(1.0, 0.0)), 0.0,
        NdcNumber(Complex128(1.0, 0.0), Complex128(0.0, 0.0))),
    PowRealSpecialTest(NdcNumber(Complex128(Float64.Inf, 0.0), Complex128(1.0, 0.0)), negZero,
        NdcNumber(Complex128(1.0, 0.0), Complex128(0.0, 0.0))),

    // PowReal(1+xϵ, y) = (1+xyϵ) for any y
    PowRealSpecialTest(NdcNumber(Complex128(1.0, 0.0), Complex128(0.0, 0.0)), 0.0,
        NdcNumber(Complex128(1.0, 0.0), Complex128(0.0, 0.0))),
    PowRealSpecialTest(NdcNumber(Complex128(1.0, 0.0), Complex128(0.0, 0.0)), 1.0,
        NdcNumber(Complex128(1.0, 0.0), Complex128(0.0, 0.0))),
    PowRealSpecialTest(NdcNumber(Complex128(1.0, 0.0), Complex128(0.0, 0.0)), 2.0,
        NdcNumber(Complex128(1.0, 0.0), Complex128(0.0, 0.0))),
    PowRealSpecialTest(NdcNumber(Complex128(1.0, 0.0), Complex128(0.0, 0.0)), 3.0,
        NdcNumber(Complex128(1.0, 0.0), Complex128(0.0, 0.0))),
    PowRealSpecialTest(NdcNumber(Complex128(1.0, 0.0), Complex128(1.0, 0.0)), 0.0,
        NdcNumber(Complex128(1.0, 0.0), Complex128(0.0, 0.0))),
    PowRealSpecialTest(NdcNumber(Complex128(1.0, 0.0), Complex128(1.0, 0.0)), 1.0,
        NdcNumber(Complex128(1.0, 0.0), Complex128(1.0, 0.0))),
    PowRealSpecialTest(NdcNumber(Complex128(1.0, 0.0), Complex128(1.0, 0.0)), 2.0,
        NdcNumber(Complex128(1.0, 0.0), Complex128(2.0, 0.0))),
    PowRealSpecialTest(NdcNumber(Complex128(1.0, 0.0), Complex128(1.0, 0.0)), 3.0,
        NdcNumber(Complex128(1.0, 0.0), Complex128(3.0, 0.0))),
    PowRealSpecialTest(NdcNumber(Complex128(1.0, 0.0), Complex128(2.0, 0.0)), 0.0,
        NdcNumber(Complex128(1.0, 0.0), Complex128(0.0, 0.0))),
    PowRealSpecialTest(NdcNumber(Complex128(1.0, 0.0), Complex128(2.0, 0.0)), 1.0,
        NdcNumber(Complex128(1.0, 0.0), Complex128(2.0, 0.0))),
    PowRealSpecialTest(NdcNumber(Complex128(1.0, 0.0), Complex128(2.0, 0.0)), 2.0,
        NdcNumber(Complex128(1.0, 0.0), Complex128(4.0, 0.0))),
    PowRealSpecialTest(NdcNumber(Complex128(1.0, 0.0), Complex128(2.0, 0.0)), 3.0,
        NdcNumber(Complex128(1.0, 0.0), Complex128(6.0, 0.0))),

    // Pow(Inf, y) = +Inf+NaNϵ for y > 0
    PowRealSpecialTest(NdcNumber(Complex128.Inf, Complex128(0.0, 0.0)), 0.5, NdcNumber(Complex128.Inf, Complex128.NaN)),
    PowRealSpecialTest(NdcNumber(Complex128.Inf, Complex128(0.0, 0.0)), 1.0, NdcNumber(Complex128.Inf, Complex128.NaN)),
    PowRealSpecialTest(NdcNumber(Complex128.Inf, Complex128(0.0, 0.0)), 1.1, NdcNumber(Complex128.Inf, Complex128.NaN)),
    PowRealSpecialTest(NdcNumber(Complex128.Inf, Complex128(0.0, 0.0)), 2.0, NdcNumber(Complex128.Inf, Complex128.NaN)),
    // Pow(Inf, y) = +0+NaNϵ for y < 0
    PowRealSpecialTest(NdcNumber(Complex128.Inf, Complex128(0.0, 0.0)), -0.5, NdcNumber(Complex128(0.0, 0.0), Complex128.NaN)),
    PowRealSpecialTest(NdcNumber(Complex128.Inf, Complex128(0.0, 0.0)), -1.0, NdcNumber(Complex128(0.0, 0.0), Complex128.NaN)),
    PowRealSpecialTest(NdcNumber(Complex128.Inf, Complex128(0.0, 0.0)), -1.1, NdcNumber(Complex128(0.0, 0.0), Complex128.NaN)),
    PowRealSpecialTest(NdcNumber(Complex128.Inf, Complex128(0.0, 0.0)), -2.0, NdcNumber(Complex128(0.0, 0.0), Complex128.NaN)),

    // PowReal(x, 1) = x for any x
    PowRealSpecialTest(NdcNumber(Complex128(0.0, 0.0), Complex128(0.0, 0.0)), 1.0,
        NdcNumber(Complex128(0.0, 0.0), Complex128(0.0, 0.0))),
    PowRealSpecialTest(NdcNumber(negZeroCmplx, Complex128(0.0, 0.0)), 1.0, NdcNumber(negZeroCmplx, Complex128(0.0, 0.0))),
    PowRealSpecialTest(NdcNumber(Complex128(0.0, 0.0), Complex128(1.0, 0.0)), 1.0,
        NdcNumber(Complex128(0.0, 0.0), Complex128(1.0, 0.0))),
    PowRealSpecialTest(NdcNumber(negZeroCmplx, Complex128(1.0, 0.0)), 1.0, NdcNumber(negZeroCmplx, Complex128(1.0, 0.0))),
    PowRealSpecialTest(NdcNumber(Complex128.NaN, Complex128(0.0, 0.0)), 1.0, NdcNumber(Complex128.NaN, Complex128(0.0, 0.0))),
    PowRealSpecialTest(NdcNumber(Complex128.NaN, Complex128(1.0, 0.0)), 1.0, NdcNumber(Complex128.NaN, Complex128(1.0, 0.0))),
    PowRealSpecialTest(NdcNumber(Complex128.NaN, Complex128(2.0, 0.0)), 1.0, NdcNumber(Complex128.NaN, Complex128(2.0, 0.0))),

    // PowReal(NaN+xϵ, y) = NaN+NaNϵ
    PowRealSpecialTest(NdcNumber(Complex128.NaN, Complex128(0.0, 0.0)), 2.0, NdcNumber(Complex128.NaN, Complex128.NaN)),
    PowRealSpecialTest(NdcNumber(Complex128.NaN, Complex128(0.0, 0.0)), 3.0, NdcNumber(Complex128.NaN, Complex128.NaN)),
    PowRealSpecialTest(NdcNumber(Complex128.NaN, Complex128(1.0, 0.0)), 2.0, NdcNumber(Complex128.NaN, Complex128.NaN)),
    PowRealSpecialTest(NdcNumber(Complex128.NaN, Complex128(1.0, 0.0)), 3.0, NdcNumber(Complex128.NaN, Complex128.NaN)),
    PowRealSpecialTest(NdcNumber(Complex128.NaN, Complex128(2.0, 0.0)), 2.0, NdcNumber(Complex128.NaN, Complex128.NaN)),
    PowRealSpecialTest(NdcNumber(Complex128.NaN, Complex128(2.0, 0.0)), 3.0, NdcNumber(Complex128.NaN, Complex128.NaN)),

    // PowReal(x, NaN) = NaN+NaNϵ
    PowRealSpecialTest(NdcNumber(Complex128(0.0, 0.0), Complex128(0.0, 0.0)), Float64.NaN,
        NdcNumber(Complex128.NaN, Complex128.NaN)),
    PowRealSpecialTest(NdcNumber(Complex128(2.0, 0.0), Complex128(0.0, 0.0)), Float64.NaN,
        NdcNumber(Complex128.NaN, Complex128.NaN)),
    PowRealSpecialTest(NdcNumber(Complex128(3.0, 0.0), Complex128(0.0, 0.0)), Float64.NaN,
        NdcNumber(Complex128.NaN, Complex128.NaN)),
    PowRealSpecialTest(NdcNumber(Complex128(0.0, 0.0), Complex128(1.0, 0.0)), Float64.NaN,
        NdcNumber(Complex128.NaN, Complex128.NaN)),
    PowRealSpecialTest(NdcNumber(Complex128(2.0, 0.0), Complex128(1.0, 0.0)), Float64.NaN,
        NdcNumber(Complex128.NaN, Complex128.NaN)),
    PowRealSpecialTest(NdcNumber(Complex128(3.0, 0.0), Complex128(1.0, 0.0)), Float64.NaN,
        NdcNumber(Complex128.NaN, Complex128.NaN)),
    PowRealSpecialTest(NdcNumber(Complex128(0.0, 0.0), Complex128(2.0, 0.0)), Float64.NaN,
        NdcNumber(Complex128.NaN, Complex128.NaN)),
    PowRealSpecialTest(NdcNumber(Complex128(2.0, 0.0), Complex128(2.0, 0.0)), Float64.NaN,
        NdcNumber(Complex128.NaN, Complex128.NaN)),
    PowRealSpecialTest(NdcNumber(Complex128(3.0, 0.0), Complex128(2.0, 0.0)), Float64.NaN,
        NdcNumber(Complex128.NaN, Complex128.NaN)),

    // Pow(-1, ±Inf) = 1
    PowRealSpecialTest(NdcNumber(Complex128(-1.0, 0.0), Complex128(0.0, 0.0)), -Float64.Inf,
        NdcNumber(Complex128(1.0, 0.0), Complex128.NaN)),
    PowRealSpecialTest(NdcNumber(Complex128(-1.0, 0.0), Complex128(0.0, 0.0)), Float64.Inf,
        NdcNumber(Complex128(1.0, 0.0), Complex128.NaN)),

    // The following tests described for cmplx.Pow ar enot valid for this type and
    // are handled by the special cases Pow(0+xϵ, y) above.
    // Pow(±0, y) = ±Inf for y an odd integer < 0
    // Pow(±0, -Inf) = +Inf
    // Pow(±0, +Inf) = +0
    // Pow(±0, y) = +Inf for finite y < 0 and not an odd integer
    // Pow(±0, y) = ±0 for y an odd integer > 0
    // Pow(±0, y) = +0 for finite y > 0 and not an odd integer

    // PowReal(x+0ϵ, +Inf) = +Inf+NaNϵ for |x| > 1
    PowRealSpecialTest(NdcNumber(Complex128(2.0, 0.0), Complex128(0.0, 0.0)), Float64.Inf,
        NdcNumber(Complex128.Inf, Complex128.NaN)),
    PowRealSpecialTest(NdcNumber(Complex128(3.0, 0.0), Complex128(0.0, 0.0)), Float64.Inf,
        NdcNumber(Complex128.Inf, Complex128.NaN)),

    // PowReal(x+yϵ, +Inf) = +Inf for |x| > 1
    PowRealSpecialTest(NdcNumber(Complex128(2.0, 0.0), Complex128(1.0, 0.0)), Float64.Inf,
        NdcNumber(Complex128.Inf, Complex128.Inf)),
    PowRealSpecialTest(NdcNumber(Complex128(3.0, 0.0), Complex128(1.0, 0.0)), Float64.Inf,
        NdcNumber(Complex128.Inf, Complex128.Inf)),
    PowRealSpecialTest(NdcNumber(Complex128(2.0, 0.0), Complex128(2.0, 0.0)), Float64.Inf,
        NdcNumber(Complex128.Inf, Complex128.Inf)),
    PowRealSpecialTest(NdcNumber(Complex128(3.0, 0.0), Complex128(2.0, 0.0)), Float64.Inf,
        NdcNumber(Complex128.Inf, Complex128.Inf)),

    // PowReal(x, -Inf) = +0+NaNϵ for |x| > 1
    PowRealSpecialTest(NdcNumber(Complex128(2.0, 0.0), Complex128(0.0, 0.0)), -Float64.Inf,
        NdcNumber(Complex128(0.0, 0.0), Complex128.NaN)),
    PowRealSpecialTest(NdcNumber(Complex128(3.0, 0.0), Complex128(0.0, 0.0)), -Float64.Inf,
        NdcNumber(Complex128(0.0, 0.0), Complex128.NaN)),
    PowRealSpecialTest(NdcNumber(Complex128(2.0, 0.0), Complex128(1.0, 0.0)), -Float64.Inf,
        NdcNumber(Complex128(0.0, 0.0), Complex128.NaN)),
    PowRealSpecialTest(NdcNumber(Complex128(3.0, 0.0), Complex128(1.0, 0.0)), -Float64.Inf,
        NdcNumber(Complex128(0.0, 0.0), Complex128.NaN)),
    PowRealSpecialTest(NdcNumber(Complex128(2.0, 0.0), Complex128(2.0, 0.0)), -Float64.Inf,
        NdcNumber(Complex128(0.0, 0.0), Complex128.NaN)),
    PowRealSpecialTest(NdcNumber(Complex128(3.0, 0.0), Complex128(2.0, 0.0)), -Float64.Inf,
        NdcNumber(Complex128(0.0, 0.0), Complex128.NaN)),

    // PowReal(x+yϵ, +Inf) = +0+NaNϵ for |x| < 1
    PowRealSpecialTest(NdcNumber(Complex128(0.1, 0.0), Complex128(0.0, 0.0)), Float64.Inf,
        NdcNumber(Complex128(0.0, 0.0), Complex128.NaN)),
    PowRealSpecialTest(NdcNumber(Complex128(0.1, 0.0), Complex128(0.1, 0.0)), Float64.Inf,
        NdcNumber(Complex128(0.0, 0.0), Complex128.NaN)),
    PowRealSpecialTest(NdcNumber(Complex128(0.2, 0.0), Complex128(0.2, 0.0)), Float64.Inf,
        NdcNumber(Complex128(0.0, 0.0), Complex128.NaN)),
    PowRealSpecialTest(NdcNumber(Complex128(0.5, 0.0), Complex128(0.5, 0.0)), Float64.Inf,
        NdcNumber(Complex128(0.0, 0.0), Complex128.NaN)),

    // PowReal(x+0ϵ, -Inf) = +Inf+NaNϵ for |x| < 1
    PowRealSpecialTest(NdcNumber(Complex128(0.1, 0.0), Complex128(0.0, 0.0)), -Float64.Inf,
        NdcNumber(Complex128.Inf, Complex128.NaN)),
    PowRealSpecialTest(NdcNumber(Complex128(0.2, 0.0), Complex128(0.0, 0.0)), -Float64.Inf,
        NdcNumber(Complex128.Inf, Complex128.NaN)),

    // PowReal(x, -Inf) = +Inf-Infϵ for |x| < 1
    PowRealSpecialTest(NdcNumber(Complex128(0.1, 0.0), Complex128(0.1, 0.0)), -Float64.Inf,
        NdcNumber(Complex128.Inf, Complex128.Inf)),
    PowRealSpecialTest(NdcNumber(Complex128(0.2, 0.0), Complex128(0.1, 0.0)), -Float64.Inf,
        NdcNumber(Complex128.Inf, Complex128.Inf)),
    PowRealSpecialTest(NdcNumber(Complex128(0.1, 0.0), Complex128(0.2, 0.0)), -Float64.Inf,
        NdcNumber(Complex128.Inf, Complex128.Inf)),
    PowRealSpecialTest(NdcNumber(Complex128(0.2, 0.0), Complex128(0.2, 0.0)), -Float64.Inf,
        NdcNumber(Complex128.Inf, Complex128.Inf)),
    PowRealSpecialTest(NdcNumber(Complex128(0.1, 0.0), Complex128(1.0, 0.0)), -Float64.Inf,
        NdcNumber(Complex128.Inf, Complex128.Inf)),
    PowRealSpecialTest(NdcNumber(Complex128(0.2, 0.0), Complex128(1.0, 0.0)), -Float64.Inf,
        NdcNumber(Complex128.Inf, Complex128.Inf)),
    PowRealSpecialTest(NdcNumber(Complex128(0.1, 0.0), Complex128(2.0, 0.0)), -Float64.Inf,
        NdcNumber(Complex128.Inf, Complex128.Inf)),
    PowRealSpecialTest(NdcNumber(Complex128(0.2, 0.0), Complex128(2.0, 0.0)), -Float64.Inf,
        NdcNumber(Complex128.Inf, Complex128.Inf))
]

@Test
public func testPowRealSpecial(): Unit {
    let tol = 1e-15
    for (test in powRealSpecialTests) {
        let got = ndcPowReal(test.d, test.p)
        @Expect(sameDual(got, test.want, tol))
    }
}

private class PowRealTest {
    let d: NdcNumber
    let p: Float64
    public init(d: NdcNumber, p: Float64) {
        this.d = d
        this.p = p
    }
}

private let powRealTests: Array<PowRealTest> = [
    PowRealTest(NdcNumber(Complex128(0.1, 0.0), Complex128(2.0, 2.0)), 0.2),
    PowRealTest(NdcNumber(Complex128(0.1, 0.0), Complex128(2.0, 2.0)), 5.0),
    PowRealTest(NdcNumber(Complex128(0.01, 0.0), Complex128(2.0, 2.0)), 0.2),
    PowRealTest(NdcNumber(Complex128(0.01, 0.0), Complex128(2.0, 2.0)), 5.0),
    PowRealTest(NdcNumber(Complex128(0.001, 0.0), Complex128(2.0, 2.0)), 0.2),
    PowRealTest(NdcNumber(Complex128(0.001, 0.0), Complex128(2.0, 2.0)), 5.0),
    PowRealTest(NdcNumber(Complex128(0.0001, 0.0), Complex128(2.0, 2.0)), 0.2),
    PowRealTest(NdcNumber(Complex128(0.0001, 0.0), Complex128(2.0, 2.0)), 5.0),
    PowRealTest(NdcNumber(Complex128(2.0, 0.0), Complex128(0.0, 0.0)), 0.5),
    PowRealTest(NdcNumber(Complex128(2.0, 0.0), Complex128(0.0, 0.0)), 2.0),
    PowRealTest(NdcNumber(Complex128(4.0, 0.0), Complex128(0.0, 0.0)), 0.5),
    PowRealTest(NdcNumber(Complex128(4.0, 0.0), Complex128(0.0, 0.0)), 2.0),
    PowRealTest(NdcNumber(Complex128(8.0, 0.0), Complex128(0.0, 0.0)), 1.0 / 3.0),
    PowRealTest(NdcNumber(Complex128(8.0, 0.0), Complex128(0.0, 0.0)), 3.0)
]

@Test
public func testPowReal(): Unit {
    let tol = 1e-14
    for (test in powRealTests) {
        var got = ndcPowReal(ndcPowReal(test.d, test.p), 1.0 / test.p)
        @Expect(sameDual(got, test.d, tol))
        if (test.p != floor(test.p)) {
            continue
        }
        let root = ndcPowReal(test.d, 1.0 / test.p)
        got = NdcNumber(Complex128(1.0, 0.0), Complex128(0.0, 0.0))
        for (_ in 0..Int64(test.p)) {
            got = ndcMul(got, root)
        }
        @Expect(sameDual(got, test.d, tol))
    }
}

private func sameDual(a: NdcNumber, b: NdcNumber, tol: Float64): Bool {
    return same(a.real, b.real, tol) && same(a.dual, b.dual, tol)
}

private func same(a: Complex128, b: Complex128, tol: Float64): Bool {
    return ((a.real.isNaN() && b.real.isNaN()) || equalWithinAbsOrRel(a.real, b.real, tol, tol)) && ((a.imag.isNaN() &&
        b.imag.isNaN()) || equalWithinAbsOrRel(a.imag, b.imag, tol, tol))
}
