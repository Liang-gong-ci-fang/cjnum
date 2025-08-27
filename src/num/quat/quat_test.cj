package cjnum.num.quat

import cjnum.util.*
import cjnum.floats.scalar.*
import std.unittest.*
import std.unittest.testmacro.*

/*
因为仓颉语言不允许出现文件之间循环调用，
所以这里把inf和nan变量放到quat_test.cj文件下
 */
public let inf = Float64.Inf
public let nan = Float64.NaN

private class ArithTest {
    let x: NqNumber
    let y: NqNumber
    let f: Float64
    let wantAdd: NqNumber
    let wantSub: NqNumber
    let wantMul: NqNumber
    let wantScale: NqNumber
    public init(x: NqNumber, y: NqNumber, f: Float64, wantAdd: NqNumber, wantSub: NqNumber, wantMul: NqNumber, wantScale: NqNumber) {
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
        NqNumber(1.0, 1.0, 1.0, 1.0),
        NqNumber(1.0, 1.0, 1.0, 1.0),
        2.0,
        NqNumber(2.0, 2.0, 2.0, 2.0),
        NqNumber(0.0, 0.0, 0.0, 0.0),
        NqNumber(-2.0, 2.0, 2.0, 2.0),
        NqNumber(2.0, 2.0, 2.0, 2.0),
    ),
    ArithTest(
        NqNumber(1.0, 1.0, 1.0, 1.0),
        NqNumber(2.0, -1.0, 1.0, -1.0),
        -2.0,
        NqNumber(3.0, 0.0, 2.0, 0.0),
        NqNumber(-1.0, 2.0, 0.0, 2.0),
        NqNumber(3.0, -1.0, 3.0, 3.0),
        NqNumber(-2.0, -2.0, -2.0, -2.0),
    ),
    ArithTest(
        NqNumber(1.0, 2.0, 3.0, 4.0),
        NqNumber(4.0, -3.0, 2.0, -1.0),
        2.0,
        NqNumber(5.0, -1.0, 5.0, 3.0),
        NqNumber(-3.0, 5.0, 1.0, 5.0),
        NqNumber(8.0, -6.0, 4.0, 28.0),
        NqNumber(2.0, 4.0, 6.0, 8.0),
    ),
    ArithTest(
        NqNumber(1.0, 2.0, 3.0, 4.0),
        NqNumber(-4.0, 3.0, -2.0, 1.0),
        -2.0,
        NqNumber(-3.0, 5.0, 1.0, 5.0),
        NqNumber(5.0, -1.0, 5.0, 3.0),
        NqNumber(-8.0, 6.0, -4.0, -28.0),
        NqNumber(-2.0, -4.0, -6.0, -8.0),
    ),
    ArithTest(
        NqNumber(-4.0, 3.0, -2.0, 1.0),
        NqNumber(1.0, 2.0, 3.0, 4.0),
        0.5,
        NqNumber(-3.0, 5.0, 1.0, 5.0),
        NqNumber(-5.0, 1.0, -5.0, -3.0),
        NqNumber(-8.0, -16.0, -24.0, -2.0),
        NqNumber(-2.0, 1.5, -1.0, 0.5),
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
    let want: NqNumber
    let wantErr: ParseError
    public init(s: String, want: NqNumber, wantErr: ParseError) {
        this.s = s
        this.want = want
        this.wantErr = wantErr
    }
}

private let parseTests: Array<ParseTest> = [
    // Simple error states:
    ParseTest("", NqNumber(0.0, 0.0, 0.0, 0.0), ParseError("", -1, r'0')),
    ParseTest("()", NqNumber(0.0, 0.0, 0.0, 0.0), ParseError("()", -1, r'0')),
    ParseTest("(1", NqNumber(0.0, 0.0, 0.0, 0.0), ParseError("(1", -1, r'0')),
    ParseTest("1)", NqNumber(0.0, 0.0, 0.0, 0.0), ParseError("1)", -1, r'0')),

    // Ambiguous parse error states:
    ParseTest("1+2i+3i", NqNumber(0.0, 0.0, 0.0, 0.0), ParseError("1+2i+3i", -1, r'0')),
    ParseTest("1+2i3j", NqNumber(0.0, 0.0, 0.0, 0.0), ParseError("1+2i3j", -1, r'0')),
    ParseTest("1e-4i-4k+10.3e6j+", NqNumber(0.0, 0.0, 0.0, 0.0), ParseError("1e-4i-4k+10.3e6j+", -1, r'0')),
    ParseTest("1e-4i-4k+10.3e6j-", NqNumber(0.0, 0.0, 0.0, 0.0), ParseError("1e-4i-4k+10.3e6j-", -1, r'0')),

    // Valid input:
    ParseTest("1+4i", NqNumber(1.0, 4.0, 0.0, 0.0), ParseError("", 0, r'0')),
    ParseTest("4i+1", NqNumber(1.0, 4.0, 0.0, 0.0), ParseError("", 0, r'0')),
    ParseTest("+1+4i", NqNumber(1.0, 4.0, 0.0, 0.0), ParseError("", 0, r'0')),
    ParseTest("+4i+1", NqNumber(1.0, 4.0, 0.0, 0.0), ParseError("", 0, r'0')),
    ParseTest("1e-4-4k+10.3e6j+1i", NqNumber(1e-4, 1.0, 10.3e6, -4.0), ParseError("", 0, r'0')),
    ParseTest("1e-4-4k+10.3e6j+i", NqNumber(1e-4, 1.0, 10.3e6, -4.0), ParseError("", 0, r'0')),
    ParseTest("1e-4-4k+10.3e6j-i", NqNumber(1e-4, -1.0, 10.3e6, -4.0), ParseError("", 0, r'0')),
    ParseTest("1e-4i-4k+10.3e6j-1", NqNumber(-1.0, 1e-4, 10.3e6, -4.0), ParseError("", 0, r'0')),
    ParseTest("1e-4i-4k+10.3e6j+1", NqNumber(1.0, 1e-4, 10.3e6, -4.0), ParseError("", 0, r'0')),
    ParseTest("(1+4i)", NqNumber(1.0, 4.0, 0.0, 0.0), ParseError("", 0, r'0')),
    ParseTest("(4i+1)", NqNumber(1.0, 4.0, 0.0, 0.0), ParseError("", 0, r'0')),
    ParseTest("(+1+4i)", NqNumber(1.0, 4.0, 0.0, 0.0), ParseError("", 0, r'0')),
    ParseTest("(+4i+1)", NqNumber(1.0, 4.0, 0.0, 0.0), ParseError("", 0, r'0')),
    ParseTest("(1e-4-4k+10.3e6j+1i)", NqNumber(1e-4, 1.0, 10.3e6, -4.0), ParseError("", 0, r'0')),
    ParseTest("(1e-4-4k+10.3e6j+i)", NqNumber(1e-4, 1.0, 10.3e6, -4.0), ParseError("", 0, r'0')),
    ParseTest("(1e-4-4k+10.3e6j-i)", NqNumber(1e-4, -1.0, 10.3e6, -4.0), ParseError("", 0, r'0')),
    ParseTest("(1e-4i-4k+10.3e6j-1)", NqNumber(-1.0, 1e-4, 10.3e6, -4.0), ParseError("", 0, r'0')),
    ParseTest("(1e-4i-4k+10.3e6j+1)", NqNumber(1.0, 1e-4, 10.3e6, -4.0), ParseError("", 0, r'0')),
    ParseTest("NaN", nqNaN(), ParseError("", 0, r'0')),
    ParseTest("nan", nqNaN(), ParseError("", 0, r'0')),
    ParseTest("Inf", nqInf(), ParseError("", 0, r'0')),
    ParseTest("inf", nqInf(), ParseError("", 0, r'0')),
    ParseTest("(Inf+Infi)", NqNumber(Float64.Inf, Float64.Inf, 0.0, 0.0), ParseError("", 0, r'0')),
    ParseTest("(-Inf+Infi)", NqNumber(-Float64.Inf, Float64.Inf, 0.0, 0.0), ParseError("", 0, r'0')),
    ParseTest("(+Inf-Infi)", NqNumber(Float64.Inf, -Float64.Inf, 0.0, 0.0), ParseError("", 0, r'0')),
    ParseTest("(inf+infi)", NqNumber(Float64.Inf, Float64.Inf, 0.0, 0.0), ParseError("", 0, r'0')),
    ParseTest("(-inf+infi)", NqNumber(-Float64.Inf, Float64.Inf, 0.0, 0.0), ParseError("", 0, r'0')),
    ParseTest("(+inf-infi)", NqNumber(Float64.Inf, -Float64.Inf, 0.0, 0.0), ParseError("", 0, r'0')),
    ParseTest("(nan+nani)", NqNumber(Float64.NaN, Float64.NaN, 0.0, 0.0), ParseError("", 0, r'0')),
    ParseTest("(nan-nani)", NqNumber(Float64.NaN, Float64.NaN, 0.0, 0.0), ParseError("", 0, r'0')),
    ParseTest("(nan+nani+1k)", NqNumber(Float64.NaN, Float64.NaN, 0.0, 1.0), ParseError("", 0, r'0')),
    ParseTest("(nan-nani+1k)", NqNumber(Float64.NaN, Float64.NaN, 0.0, 1.0), ParseError("", 0, r'0'))
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

public func nqEqualApprox(a: NqNumber, b: NqNumber, tol: Float64): Bool {
    return equalWithinAbsOrRel(a.real, b.real, tol, tol) && equalWithinAbsOrRel(a.imag, b.imag, tol, tol) &&
        equalWithinAbsOrRel(a.jmag, b.jmag, tol, tol) && equalWithinAbsOrRel(a.kmag, b.kmag, tol, tol)
}

public func nqSameApprox(a: NqNumber, b: NqNumber, tol: Float64): Bool {
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

private func sameNumber(a: NqNumber, b: NqNumber): Bool {
    return sameFloat(a.real, b.real) && sameFloat(a.imag, b.imag) && sameFloat(a.jmag, b.jmag) && sameFloat(a.kmag,
        b.kmag)
}

private func sameFloat(a: Float64, b: Float64): Bool {
    return a == b || (a.isNaN() && b.isNaN())
}
