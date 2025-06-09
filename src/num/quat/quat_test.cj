package cjnum.num.quat

import cjnum.util.*
import cjnum.floats.scalar.*
import std.math.*
import std.unittest.*
import std.unittest.testmacro.*

/*
因为仓颉语言不允许出现文件之间循环调用，
所以这里把inf和nan变量放到quat_test.cj文件下
 */
public let inf = Float64.Inf
public let nan = Float64.NaN

private class ArithTest {
    let x: Number
    let y: Number
    let f: Float64
    let wantAdd: Number
    let wantSub: Number
    let wantMul: Number
    let wantScale: Number
    public init(x: Number, y: Number, f: Float64, wantAdd: Number, wantSub: Number, wantMul: Number, wantScale: Number) {
        this.x = x
        this.y = y
        this.f = f
        this.wantAdd = wantAdd
        this.wantSub = wantSub
        this.wantMul = wantMul
        this.wantScale = wantScale
    }
}

private let arithTests: Array<ArithTest> = [
    ArithTest(
        Number(1.0, 1.0, 1.0, 1.0),
        Number(1.0, 1.0, 1.0, 1.0),
        2.0,
        Number(2.0, 2.0, 2.0, 2.0),
        Number(0.0, 0.0, 0.0, 0.0),
        Number(-2.0, 2.0, 2.0, 2.0),
        Number(2.0, 2.0, 2.0, 2.0),
    ),
    ArithTest(
        Number(1.0, 1.0, 1.0, 1.0),
        Number(2.0, -1.0, 1.0, -1.0),
        -2.0,
        Number(3.0, 0.0, 2.0, 0.0),
        Number(-1.0, 2.0, 0.0, 2.0),
        Number(3.0, -1.0, 3.0, 3.0),
        Number(-2.0, -2.0, -2.0, -2.0),
    ),
    ArithTest(
        Number(1.0, 2.0, 3.0, 4.0),
        Number(4.0, -3.0, 2.0, -1.0),
        2.0,
        Number(5.0, -1.0, 5.0, 3.0),
        Number(-3.0, 5.0, 1.0, 5.0),
        Number(8.0, -6.0, 4.0, 28.0),
        Number(2.0, 4.0, 6.0, 8.0),
    ),
    ArithTest(
        Number(1.0, 2.0, 3.0, 4.0),
        Number(-4.0, 3.0, -2.0, 1.0),
        -2.0,
        Number(-3.0, 5.0, 1.0, 5.0),
        Number(5.0, -1.0, 5.0, 3.0),
        Number(-8.0, 6.0, -4.0, -28.0),
        Number(-2.0, -4.0, -6.0, -8.0),
    ),
    ArithTest(
        Number(-4.0, 3.0, -2.0, 1.0),
        Number(1.0, 2.0, 3.0, 4.0),
        0.5,
        Number(-3.0, 5.0, 1.0, 5.0),
        Number(-5.0, 1.0, -5.0, -3.0),
        Number(-8.0, -16.0, -24.0, -2.0),
        Number(-2.0, 1.5, -1.0, 0.5),
    )
]

@Test
public func testArithmetic(): Unit {
    for (test in arithTests) {
        let gotAdd = nqAdd(test.x, test.y)
        @Expect(gotAdd.nqEqual(test.wantAdd))
        let gotSub = nqSub(test.x, test.y)
        @Expect(gotSub.nqEqual(test.wantSub))
        let gotMul = nqMul(test.x, test.y)
        @Expect(gotMul.nqEqual(test.wantMul))
        let gotScale = nqScale(test.f, test.x)
        @Expect(gotScale.nqEqual(test.wantScale))
    }
}

// TODO var formatTests = []struct {

// TODO func TestFormat(t *testing.T) {

private class ParseTest {
    let s: String
    let want: Number
    let wantErr: ParseError
    public init(s: String, want: Number, wantErr: ParseError) {
        this.s = s
        this.want = want
        this.wantErr = wantErr
    }
}

private let parseTests: Array<ParseTest> = [
    // Simple error states:
    ParseTest("", Number(0.0, 0.0, 0.0, 0.0), ParseError("", -1, r'0')),
    ParseTest("()", Number(0.0, 0.0, 0.0, 0.0), ParseError("()", -1, r'0')),
    ParseTest("(1", Number(0.0, 0.0, 0.0, 0.0), ParseError("(1", -1, r'0')),
    ParseTest("1)", Number(0.0, 0.0, 0.0, 0.0), ParseError("1)", -1, r'0')),

    // Ambiguous parse error states:
    ParseTest("1+2i+3i", Number(0.0, 0.0, 0.0, 0.0), ParseError("1+2i+3i", -1, r'0')),
    ParseTest("1+2i3j", Number(0.0, 0.0, 0.0, 0.0), ParseError("1+2i3j", -1, r'0')),
    ParseTest("1e-4i-4k+10.3e6j+", Number(0.0, 0.0, 0.0, 0.0), ParseError("1e-4i-4k+10.3e6j+", -1, r'0')),
    ParseTest("1e-4i-4k+10.3e6j-", Number(0.0, 0.0, 0.0, 0.0), ParseError("1e-4i-4k+10.3e6j-", -1, r'0')),

    // Valid input:
    ParseTest("1+4i", Number(1.0, 4.0, 0.0, 0.0), ParseError("", 0, r'0')),
    ParseTest("4i+1", Number(1.0, 4.0, 0.0, 0.0), ParseError("", 0, r'0')),
    ParseTest("+1+4i", Number(1.0, 4.0, 0.0, 0.0), ParseError("", 0, r'0')),
    ParseTest("+4i+1", Number(1.0, 4.0, 0.0, 0.0), ParseError("", 0, r'0')),
    ParseTest("1e-4-4k+10.3e6j+1i", Number(1e-4, 1.0, 10.3e6, -4.0), ParseError("", 0, r'0')),
    ParseTest("1e-4-4k+10.3e6j+i", Number(1e-4, 1.0, 10.3e6, -4.0), ParseError("", 0, r'0')),
    ParseTest("1e-4-4k+10.3e6j-i", Number(1e-4, -1.0, 10.3e6, -4.0), ParseError("", 0, r'0')),
    ParseTest("1e-4i-4k+10.3e6j-1", Number(-1.0, 1e-4, 10.3e6, -4.0), ParseError("", 0, r'0')),
    ParseTest("1e-4i-4k+10.3e6j+1", Number(1.0, 1e-4, 10.3e6, -4.0), ParseError("", 0, r'0')),
    ParseTest("(1+4i)", Number(1.0, 4.0, 0.0, 0.0), ParseError("", 0, r'0')),
    ParseTest("(4i+1)", Number(1.0, 4.0, 0.0, 0.0), ParseError("", 0, r'0')),
    ParseTest("(+1+4i)", Number(1.0, 4.0, 0.0, 0.0), ParseError("", 0, r'0')),
    ParseTest("(+4i+1)", Number(1.0, 4.0, 0.0, 0.0), ParseError("", 0, r'0')),
    ParseTest("(1e-4-4k+10.3e6j+1i)", Number(1e-4, 1.0, 10.3e6, -4.0), ParseError("", 0, r'0')),
    ParseTest("(1e-4-4k+10.3e6j+i)", Number(1e-4, 1.0, 10.3e6, -4.0), ParseError("", 0, r'0')),
    ParseTest("(1e-4-4k+10.3e6j-i)", Number(1e-4, -1.0, 10.3e6, -4.0), ParseError("", 0, r'0')),
    ParseTest("(1e-4i-4k+10.3e6j-1)", Number(-1.0, 1e-4, 10.3e6, -4.0), ParseError("", 0, r'0')),
    ParseTest("(1e-4i-4k+10.3e6j+1)", Number(1.0, 1e-4, 10.3e6, -4.0), ParseError("", 0, r'0')),
    ParseTest("NaN", nqNaN(), ParseError("", 0, r'0')),
    ParseTest("nan", nqNaN(), ParseError("", 0, r'0')),
    ParseTest("Inf", nqInf(), ParseError("", 0, r'0')),
    ParseTest("inf", nqInf(), ParseError("", 0, r'0')),
    ParseTest("(Inf+Infi)", Number(Float64.Inf, Float64.Inf, 0.0, 0.0), ParseError("", 0, r'0')),
    ParseTest("(-Inf+Infi)", Number(-Float64.Inf, Float64.Inf, 0.0, 0.0), ParseError("", 0, r'0')),
    ParseTest("(+Inf-Infi)", Number(Float64.Inf, -Float64.Inf, 0.0, 0.0), ParseError("", 0, r'0')),
    ParseTest("(inf+infi)", Number(Float64.Inf, Float64.Inf, 0.0, 0.0), ParseError("", 0, r'0')),
    ParseTest("(-inf+infi)", Number(-Float64.Inf, Float64.Inf, 0.0, 0.0), ParseError("", 0, r'0')),
    ParseTest("(+inf-infi)", Number(Float64.Inf, -Float64.Inf, 0.0, 0.0), ParseError("", 0, r'0')),
    ParseTest("(nan+nani)", Number(Float64.NaN, Float64.NaN, 0.0, 0.0), ParseError("", 0, r'0')),
    ParseTest("(nan-nani)", Number(Float64.NaN, Float64.NaN, 0.0, 0.0), ParseError("", 0, r'0')),
    ParseTest("(nan+nani+1k)", Number(Float64.NaN, Float64.NaN, 0.0, 1.0), ParseError("", 0, r'0')),
    ParseTest("(nan-nani+1k)", Number(Float64.NaN, Float64.NaN, 0.0, 1.0), ParseError("", 0, r'0'))
]

@Test
public func testParse(): Unit {
    for (test in parseTests) {
        if (test.wantErr.state == -1) {
            let _ = @ExpectThrows[Exception](nqParse(test.s))
        } else {
            let got = nqParse(test.s)
            @Expect(sameNumber(got, test.want))
        }
    }
}

public func nqEqualApprox(a: Number, b: Number, tol: Float64): Bool {
    return equalWithinAbsOrRel(a.real, b.real, tol, tol) && equalWithinAbsOrRel(a.imag, b.imag, tol, tol) &&
        equalWithinAbsOrRel(a.jmag, b.jmag, tol, tol) && equalWithinAbsOrRel(a.kmag, b.kmag, tol, tol)
}

public func nqSameApprox(a: Number, b: Number, tol: Float64): Bool {
    if (a.real == 0.0 && b.real == 0.0) {
        return signbit(a.real) == signbit(b.real)
    } else if (a.imag == 0.0 && b.imag == 0.0) {
        return signbit(a.imag) == signbit(b.imag)
    } else if (a.jmag == 0.0 && b.jmag == 0.0) {
        return signbit(a.jmag) == signbit(b.jmag)
    } else if (a.kmag == 0.0 && b.kmag == 0.0) {
        return signbit(a.kmag) == signbit(b.kmag)
    }
    return (sameFloat(a.real, b.real) || equalWithinAbsOrRel(a.real, b.real, tol, tol)) && (sameFloat(a.imag, b.imag) ||
        equalWithinAbsOrRel(a.imag, b.imag, tol, tol)) && (sameFloat(a.jmag, b.jmag) || equalWithinAbsOrRel(a.jmag,
        b.jmag, tol, tol)) && (sameFloat(a.kmag, b.kmag) || equalWithinAbsOrRel(a.kmag, b.kmag, tol, tol))
}

private func sameNumber(a: Number, b: Number): Bool {
    return sameFloat(a.real, b.real) && sameFloat(a.imag, b.imag) && sameFloat(a.jmag, b.jmag) && sameFloat(a.kmag,
        b.kmag)
}

private func sameFloat(a: Float64, b: Float64): Bool {
    return a == b || (a.isNaN() && b.isNaN())
}
