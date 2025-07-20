package cjnum.num.quat

import cjnum.util.*
import cjnum.complex.*
import std.math.*
import std.unittest.*
import std.unittest.testmacro.*

private class SinTest {
    let q: NqNumber
    let want: NqNumber
    public init(q: NqNumber, want: NqNumber) {
        this.q = q
        this.want = want
    }
}

private let sinTests: Array<SinTest> = [
    SinTest(NqNumber(0.0, 0.0, 0.0, 0.0), NqNumber(0.0, 0.0, 0.0, 0.0)),
    SinTest(NqNumber(Float64.getPI() / 2.0, 0.0, 0.0, 0.0), NqNumber(1.0, 0.0, 0.0, 0.0)),
    SinTest(
        NqNumber(0.0, Float64.getPI() / 2.0, 0.0, 0.0),
        NqNumber(0.0, c128Sin(Complex128(0.0, Float64.getPI() / 2.0)).imag, 0.0, 0.0)
    ),
    SinTest(
        NqNumber(0.0, 0.0, Float64.getPI() / 2.0, 0.0),
        NqNumber(0.0, 0.0, c128Sin(Complex128(0.0, Float64.getPI() / 2.0)).imag, 0.0)
    ),
    SinTest(
        NqNumber(0.0, 0.0, 0.0, Float64.getPI() / 2.0),
        NqNumber(0.0, 0.0, 0.0, c128Sin(Complex128(0.0, Float64.getPI() / 2.0)).imag)
    ),

    // Exercises from Real Quaternionic Calculus Handbook doi:10.1007/978-3-0348-0622-0
    // Ex 6.159 (a) and (b).
    SinTest(
        NqNumber(1.0, 1.0, 1.0, 1.0),
        NqNumber(sin(1.0) * cosh(sqrt(3.0)), cos(1.0) * sinh(sqrt(3.0)) / sqrt(3.0), cos(1.0) * sinh(sqrt(3.0)) / sqrt(
            3.0), cos(1.0) * sinh(sqrt(3.0)) / sqrt(3.0))
    ),
    SinTest(
        NqNumber(0.0, -2.0, 1.0, 0.0),
        NqNumber(0.0, -2.0 * sinh(sqrt(5.0)) / sqrt(5.0), sinh(sqrt(5.0)) / sqrt(5.0), 0.0)
    )
]

@Test
public func testSin(): Unit {
    let tol = 1e-14
    for (test in sinTests) {
        let got = nqSin(test.q)
        @Expect(nqEqualApprox(got, test.want, tol))
    }
}

private class SinhTest {
    let q: NqNumber
    let want: NqNumber
    public init(q: NqNumber, want: NqNumber) {
        this.q = q
        this.want = want
    }
}

private let sinhTests: Array<SinhTest> = [
    SinhTest(NqNumber(0.0, 0.0, 0.0, 0.0), NqNumber(0.0, 0.0, 0.0, 0.0)),
    SinhTest(NqNumber(Float64.getPI() / 2.0, 0.0, 0.0, 0.0), NqNumber(sinh(Float64.getPI() / 2.0), 0.0, 0.0, 0.0)),
    SinhTest(
        NqNumber(0.0, Float64.getPI() / 2.0, 0.0, 0.0),
        NqNumber(0.0, c128Sinh(Complex128(0.0, Float64.getPI() / 2.0)).imag, 0.0, 0.0)
    ),
    SinhTest(
        NqNumber(0.0, 0.0, Float64.getPI() / 2.0, 0.0),
        NqNumber(0.0, 0.0, c128Sinh(Complex128(0.0, Float64.getPI() / 2.0)).imag, 0.0)
    ),
    SinhTest(
        NqNumber(0.0, 0.0, 0.0, Float64.getPI() / 2.0),
        NqNumber(0.0, 0.0, 0.0, c128Sinh(Complex128(0.0, Float64.getPI() / 2.0)).imag)
    ),
    SinhTest(
        NqNumber(1.0, -1.0, -1.0, 0.0),
        nqScale(
            0.5,
            nqSub(
                nqExp(NqNumber(1.0, -1.0, -1.0, 0.0)),
                nqExp(nqScale(-1.0, NqNumber(1.0, -1.0, -1.0, 0.0)))
            )
        )
    ),
    SinhTest(
        NqNumber(1.0, 1.0, 1.0, 1.0),
        nqScale(
            0.5,
            nqSub(
                nqExp(NqNumber(1.0, 1.0, 1.0, 1.0)),
                nqExp(nqScale(-1.0, NqNumber(1.0, 1.0, 1.0, 1.0)))
            )
        )
    ),
    SinhTest(
        nqAsinh(NqNumber(1.0, 1.0, 1.0, 1.0)),
        NqNumber(1.0, 1.0, 1.0, 1.0)
    ),
    SinhTest(
        nqAsinh(NqNumber(1.0, 1.0, 1.0, 1.0)),
        nqScale(
            0.5,
            nqSub(
                nqExp(nqAsinh(NqNumber(1.0, 1.0, 1.0, 1.0))),
                nqExp(nqScale(-1.0, nqAsinh(NqNumber(1.0, 1.0, 1.0, 1.0))))
            )
        )
    ),
    SinhTest(
        NqNumber(Float64.Inf, 0.0, 0.0, 0.0),
        NqNumber(Float64.Inf, 0.0, 0.0, 0.0)
    ),
    SinhTest(
        NqNumber(Float64.Inf, Float64.getPI() / 2.0, 0.0, 0.0),
        NqNumber(Float64.Inf, Float64.Inf, 0.0, 0.0)
    ),
    SinhTest(
        NqNumber(Float64.Inf, Float64.getPI(), 0.0, 0.0),
        NqNumber(-Float64.Inf, Float64.Inf, 0.0, 0.0)
    ),
    SinhTest(
        NqNumber(Float64.Inf, 3.0 * Float64.getPI() / 2.0, 0.0, 0.0),
        NqNumber(-Float64.Inf, -Float64.Inf, 0.0, 0.0)
    ),
    SinhTest(
        NqNumber(Float64.Inf, 2.0 * Float64.getPI(), 0.0, 0.0),
        NqNumber(Float64.Inf, -Float64.Inf, 0.0, 0.0)
    )
]

@Test
public func testSinh(): Unit {
    let tol = 1e-14
    for (test in sinhTests) {
        let got = nqSinh(test.q)
        @Expect(nqSameApprox(got, test.want, tol))
    }
}

private class CosTest {
    let q: NqNumber
    let want: NqNumber
    public init(q: NqNumber, want: NqNumber) {
        this.q = q
        this.want = want
    }
}

private let cosTests: Array<CosTest> = [
    CosTest(NqNumber(0.0, 0.0, 0.0, 0.0), NqNumber(1.0, 0.0, 0.0, 0.0)),
    CosTest(
        NqNumber(Float64.getPI() / 2.0, 0.0, 0.0, 0.0),
        NqNumber(0.0, 0.0, 0.0, 0.0)
    ),
    CosTest(
        NqNumber(0.0, Float64.getPI() / 2.0, 0.0, 0.0),
        NqNumber(c128Cos(Complex128(0.0, Float64.getPI() / 2.0)).real, 0.0, 0.0, 0.0)
    ),
    CosTest(
        NqNumber(0.0, 0.0, Float64.getPI() / 2.0, 0.0),
        NqNumber((c128Cos(Complex128(0.0, Float64.getPI() / 2.0))).real, 0.0, 0.0, 0.0)
    ),
    CosTest(
        NqNumber(0.0, 0.0, 0.0, Float64.getPI() / 2.0),
        NqNumber((c128Cos(Complex128(0.0, Float64.getPI() / 2.0))).real, 0.0, 0.0, 0.0)
    ),
    CosTest(
        NqNumber(1.0, 1.0, 1.0, 1.0),
        NqNumber(
            cos(1.0) * cosh(sqrt(3.0)),
            -sin(1.0) * sinh(sqrt(3.0)) / sqrt(3.0),
            -sin(1.0) * sinh(sqrt(3.0)) / sqrt(3.0),
            -sin(1.0) * sinh(sqrt(3.0)) / sqrt(3.0)
        )
    )
]

@Test
public func testCos(): Unit {
    let tol = 1e-14
    for (test in cosTests) {
        let got = nqCos(test.q)
        @Expect(nqEqualApprox(got, test.want, tol))
    }
}

private class CoshTest {
    let q: NqNumber
    let want: NqNumber
    public init(q: NqNumber, want: NqNumber) {
        this.q = q
        this.want = want
    }
}

private let coshTests: Array<CoshTest> = [
    CoshTest(NqNumber(0.0, 0.0, 0.0, 0.0), NqNumber(1.0, 0.0, 0.0, 0.0)),
    CoshTest(
        NqNumber(Float64.getPI() / 2.0, 0.0, 0.0, 0.0),
        NqNumber(cosh(Float64.getPI() / 2.0), 0.0, 0.0, 0.0)
    ),
    CoshTest(
        NqNumber(0.0, Float64.getPI() / 2.0, 0.0, 0.0),
        NqNumber(0.0, c128Cosh(Complex128(0.0, Float64.getPI() / 2.0)).imag, 0.0, 0.0)
    ),
    CoshTest(
        NqNumber(0.0, 0.0, Float64.getPI() / 2.0, 0.0),
        NqNumber(0.0, 0.0, c128Cosh(Complex128(0.0, Float64.getPI() / 2.0)).imag, 0.0)
    ),
    CoshTest(
        NqNumber(0.0, 0.0, 0.0, Float64.getPI() / 2.0),
        NqNumber(0.0, 0.0, 0.0, c128Cosh(Complex128(0.0, Float64.getPI() / 2.0)).imag)
    ),
    CoshTest(
        NqNumber(1.0, -1.0, -1.0, 0.0),
        nqScale(0.5, nqAdd(nqExp(NqNumber(1.0, -1.0, -1.0, 0.0)), nqExp(nqScale(-1.0, NqNumber(1.0, -1.0, -1.0, 0.0)))))
    ),
    CoshTest(
        NqNumber(1.0, 1.0, 1.0, 1.0),
        nqScale(0.5, nqAdd(nqExp(NqNumber(1.0, 1.0, 1.0, 1.0)), nqExp(nqScale(-1.0, NqNumber(1.0, 1.0, 1.0, 1.0)))))
    ),
    CoshTest(
        NqNumber(Float64.Inf, 0.0, 0.0, 0.0),
        NqNumber(Float64.Inf, 0.0, 0.0, 0.0)
    ),
    CoshTest(
        NqNumber(Float64.Inf, Float64.getPI() / 2.0, 0.0, 0.0),
        NqNumber(Float64.Inf, Float64.Inf, 0.0, 0.0)
    ),
    CoshTest(
        NqNumber(Float64.Inf, Float64.getPI(), 0.0, 0.0),
        NqNumber(-Float64.Inf, Float64.Inf, 0.0, 0.0)
    ),
    CoshTest(
        NqNumber(Float64.Inf, 3.0 * Float64.getPI() / 2.0, 0.0, 0.0),
        NqNumber(-Float64.Inf, -Float64.Inf, 0.0, 0.0)
    ),
    CoshTest(
        NqNumber(Float64.Inf, 2.0 * Float64.getPI(), 0.0, 0.0),
        NqNumber(Float64.Inf, -Float64.Inf, 0.0, 0.0)
    )
]

@Test
public func testCosh(): Unit {
    let tol = 1e-14
    for (test in coshTests) {
        let got = nqCosh(test.q)
        @Expect(nqSameApprox(got, test.want, tol))
    }
}

private class TanTest {
    let q: NqNumber
    let want: NqNumber
    public init(q: NqNumber, want: NqNumber) {
        this.q = q
        this.want = want
    }
}

private let tanTests: Array<TanTest> = [
    TanTest(NqNumber(0.0, 0.0, 0.0, 0.0), NqNumber(1.0, 0.0, 0.0, 0.0)),
    TanTest(
        NqNumber(Float64.getPI() / 4.0, 0.0, 0.0, 0.0),
        NqNumber(tan(Float64.getPI() / 4.0), 0.0, 0.0, 0.0)
    ),
    TanTest(
        NqNumber(0.0, Float64.getPI() / 4.0, 0.0, 0.0),
        NqNumber(0.0, c128Tan(Complex128(0.0, Float64.getPI() / 4.0)).imag, 0.0, 0.0)
    ),
    TanTest(
        NqNumber(0.0, 0.0, Float64.getPI() / 4.0, 0.0),
        NqNumber(0.0, 0.0, c128Tan(Complex128(0.0, Float64.getPI() / 4.0)).imag, 0.0)
    ),
    TanTest(
        NqNumber(0.0, 0.0, 0.0, Float64.getPI() / 4.0),
        NqNumber(0.0, 0.0, 0.0, c128Tan(Complex128(0.0, Float64.getPI() / 4.0)).imag)
    ),
    TanTest(
        NqNumber(0.0, 1.0, 0.0, 0.0),
        nqMul(
            nqSin(NqNumber(0.0, 1.0, 0.0, 0.0)),
            nqInv(nqCos(NqNumber(0.0, 1.0, 0.0, 0.0)))
        )
    ),
    TanTest(
        NqNumber(1.0, 1.0, 1.0, 1.0),
        nqMul(
            nqSin(NqNumber(1.0, 1.0, 1.0, 1.0)),
            nqInv(nqCos(NqNumber(1.0, 1.0, 1.0, 1.0)))
        )
    )
]

public func testTan(): Unit {
    let tol = 1e-14
    for (test in tanTests) {
        let got = nqTan(test.q)
        @Expect(nqEqualApprox(got, test.want, tol))
    }
}

private class TanhTest {
    let q: NqNumber
    let want: NqNumber
    public init(q: NqNumber, want: NqNumber) {
        this.q = q
        this.want = want
    }
}

private let tanhTests: Array<TanhTest> = [
    TanhTest(NqNumber(0.0, 0.0, 0.0, 0.0), NqNumber(0.0, 0.0, 0.0, 0.0)),
    TanhTest(
        NqNumber(Float64.getPI() / 4.0, 0.0, 0.0, 0.0),
        NqNumber(tanh(Float64.getPI() / 4.0), 0.0, 0.0, 0.0)
    ),
    TanhTest(
        NqNumber(0.0, Float64.getPI() / 4.0, 0.0, 0.0),
        NqNumber(0.0, c128Tanh(Complex128(0.0, Float64.getPI() / 4.0)).imag, 0.0, 0.0)
    ),
    TanhTest(
        NqNumber(0.0, 0.0, Float64.getPI() / 4.0, 0.0),
        NqNumber(0.0, 0.0, c128Tanh(Complex128(0.0, Float64.getPI() / 4.0)).imag, 0.0)
    ),
    TanhTest(
        NqNumber(0.0, 0.0, 0.0, Float64.getPI() / 4.0),
        NqNumber(0.0, 0.0, 0.0, c128Tanh(Complex128(0.0, Float64.getPI() / 4.0)).imag)
    ),
    TanhTest(
        NqNumber(0.0, 1.0, 0.0, 0.0),
        nqMul(nqSinh(NqNumber(0.0, 1.0, 0.0, 0.0)), nqInv(nqCosh(NqNumber(0.0, 1.0, 0.0, 0.0))))
    ),
    TanhTest(
        NqNumber(1.0, 1.0, 1.0, 1.0),
        nqMul(nqSinh(NqNumber(1.0, 1.0, 1.0, 1.0)), nqInv(nqCosh(NqNumber(1.0, 1.0, 1.0, 1.0))))
    ),
    TanhTest(
        NqNumber(Float64.Inf, 0.0, 0.0, 0.0),
        NqNumber(1.0, 0.0, 0.0, 0.0)
    ),
    TanhTest(
        NqNumber(Float64.Inf, Float64.getPI() / 4.0, 0.0, 0.0),
        NqNumber(1.0, 0.0 * sin(Float64.getPI() / 2.0), 0.0, 0.0)
    ),
    TanhTest(
        NqNumber(Float64.Inf, Float64.getPI() / 2.0, 0.0, 0.0),
        NqNumber(1.0, 0.0 * sin(Float64.getPI()), 0.0, 0.0)
    ),
    TanhTest(
        NqNumber(Float64.Inf, 3.0 * Float64.getPI() / 4.0, 0.0, 0.0),
        NqNumber(1.0, 0.0 * sin(3.0 * Float64.getPI() / 2.0), 0.0, 0.0)
    ),
    TanhTest(
        NqNumber(Float64.Inf, Float64.getPI(), 0.0, 0.0),
        NqNumber(1.0, 0.0 * sin(2.0 * Float64.getPI()), 0.0, 0.0)
    )
]

@Test
public func testTanh(): Unit {
    let tol = 1e-14
    for (test in tanhTests) {
        let got = nqTanh(test.q)
        @Expect(nqSameApprox(got, test.want, tol))
    }
}

private class AsinTest {
    let q: NqNumber
    let want: NqNumber
    public init(q: NqNumber, want: NqNumber) {
        this.q = q
        this.want = want
    }
}

private let asinTests: Array<AsinTest> = [
    AsinTest(NqNumber(0.0, 0.0, 0.0, 0.0), NqNumber(0.0, 0.0, 0.0, 0.0)),
    AsinTest(
        NqNumber(1.0, 0.0, 0.0, 0.0),
        NqNumber(Float64.getPI() / 2.0, 0.0, 0.0, 0.0)
    ),
    AsinTest(
        NqNumber(0.0, 1.0, 0.0, 0.0),
        NqNumber(0.0, c128Asinh(Complex128(1.0, 0.0)).real, 0.0, 0.0)
    ),
    AsinTest(
        NqNumber(0.0, 0.0, 1.0, 0.0),
        NqNumber(0.0, 0.0, c128Asinh(Complex128(1.0, 0.0)).real, 0.0)
    ),
    AsinTest(
        NqNumber(0.0, 0.0, 0.0, 1.0),
        NqNumber(0.0, 0.0, 0.0, c128Asinh(Complex128(1.0, 0.0)).real)
    ),
    AsinTest(
        nqSin(NqNumber(1.0, 1.0, 1.0, 1.0)),
        NqNumber(1.0, 1.0, 1.0, 1.0)
    )
]

@Test
public func testAsin(): Unit {
    let tol = 1e-14
    for (test in asinTests) {
        let got = nqAsin(test.q)
        @Expect(nqEqualApprox(got, test.want, tol))
    }
}

private class AsinhTest {
    let q: NqNumber
    let want: NqNumber
    public init(q: NqNumber, want: NqNumber) {
        this.q = q
        this.want = want
    }
}

private let asinhTests: Array<AsinhTest> = [
    AsinhTest(NqNumber(0.0, 0.0, 0.0, 0.0), NqNumber(0.0, 0.0, 0.0, 0.0)),
    AsinhTest(
        NqNumber(1.0, 0.0, 0.0, 0.0),
        NqNumber(asinh(1.0), 0.0, 0.0, 0.0)
    ),
    AsinhTest(
        NqNumber(0.0, 1.0, 0.0, 0.0),
        NqNumber(0.0, Float64.getPI() / 2.0, 0.0, 0.0)
    ),
    AsinhTest(
        NqNumber(0.0, 0.0, 1.0, 0.0),
        NqNumber(0.0, 0.0, Float64.getPI() / 2.0, 0.0)
    ),
    AsinhTest(
        NqNumber(0.0, 0.0, 0.0, 1.0),
        NqNumber(0.0, 0.0, 0.0, Float64.getPI() / 2.0)
    ),
    AsinhTest(
        NqNumber(1.0, 1.0, 1.0, 1.0),
        nqLog(
            nqAdd(
                NqNumber(1.0, 1.0, 1.0, 1.0),
                nqSqrt(
                    nqAdd(
                        nqMul(NqNumber(1.0, 1.0, 1.0, 1.0), NqNumber(1.0, 1.0, 1.0, 1.0)),
                        NqNumber(1.0, 0.0, 0.0, 0.0)
                    )
                )
            )
        )
    ),
    AsinhTest(
        nqSinh(NqNumber(1.0, 0.0, 0.0, 0.0)),
        NqNumber(1.0, 0.0, 0.0, 0.0)
    ),
    AsinhTest(
        nqSinh(NqNumber(0.0, 1.0, 0.0, 0.0)),
        NqNumber(0.0, 1.0, 0.0, 0.0)
    ),
    AsinhTest(
        nqSinh(NqNumber(0.0, 1.0, 1.0, 0.0)),
        NqNumber(0.0, 1.0, 1.0, 0.0)
    ),
    AsinhTest(
        nqSinh(NqNumber(1.0, 1.0, 1.0, 0.0)),
        NqNumber(1.0, 1.0, 1.0, 0.0)
    ),
    // The following fails:
    // {q: Sinh(NqNumber{1, 1, 1, 1}), want: NqNumber{1, 1, 1, 1}},
    // but this passes...
    AsinhTest(
        nqSinh(NqNumber(1.0, 1.0, 1.0, 1.0)),
        nqLog(
            nqAdd(
                nqSinh(NqNumber(1.0, 1.0, 1.0, 1.0)),
                nqSqrt(
                    nqAdd(
                        nqMul(nqSinh(NqNumber(1.0, 1.0, 1.0, 1.0)), nqSinh(NqNumber(1.0, 1.0, 1.0, 1.0))),
                        NqNumber(1.0, 0.0, 0.0, 0.0)
                    )
                )
            )
        )
    )
    // And see the Sinh tests that do the reciprocal operation.
]

@Test
public func testAsinh() {
    let tol = 1e-14
    for (test in asinhTests) {
        let got = nqAsinh(test.q)
        @Expect(nqEqualApprox(got, test.want, tol))
    }
}

private class AcosTest {
    let q: NqNumber
    let want: NqNumber
    public init(q: NqNumber, want: NqNumber) {
        this.q = q
        this.want = want
    }
}

private let acosTests: Array<AcosTest> = [
    AcosTest(NqNumber(0.0, 0.0, 0.0, 0.0), NqNumber(Float64.getPI() / 2.0, 0.0, 0.0, 0.0)),
    AcosTest(
        NqNumber(1.0, 0.0, 0.0, 0.0),
        NqNumber(0.0, 0.0, 0.0, 0.0)
    ),
    AcosTest(
        NqNumber(0.0, 1.0, 0.0, 0.0),
        NqNumber(
            c128Acos(Complex128(0.0, 1.0)).real,
            c128Acos(Complex128(0.0, 1.0)).imag,
            0.0,
            0.0
        )
    ),
    AcosTest(
        NqNumber(0.0, 0.0, 1.0, 0.0),
        NqNumber(
            c128Acos(Complex128(0.0, 1.0)).real,
            0.0,
            c128Acos(Complex128(0.0, 1.0)).imag,
            0.0
        )
    ),
    AcosTest(
        NqNumber(0.0, 0.0, 0.0, 1.0),
        NqNumber(
            c128Acos(Complex128(0.0, 1.0)).real,
            0.0,
            0.0,
            c128Acos(Complex128(0.0, 1.0)).imag
        )
    ),
    AcosTest(
        nqCos(NqNumber(1.0, 1.0, 1.0, 1.0)),
        NqNumber(1.0, 1.0, 1.0, 1.0)
    )
]

@Test
public func testAcos(): Unit {
    let tol = 1e-14
    for (test in acosTests) {
        let got = nqAcos(test.q)
        @Expect(nqEqualApprox(got, test.want, tol))
    }
}

private class AcoshTest {
    let q: NqNumber
    let want: NqNumber
    public init(q: NqNumber, want: NqNumber) {
        this.q = q
        this.want = want
    }
}

private let acoshTests: Array<AcoshTest> = [
    AcoshTest(NqNumber(0.0, 0.0, 0.0, 0.0), NqNumber(Float64.getPI() / 2.0, 0.0, 0.0, 0.0)),
    AcoshTest(
        NqNumber(1.0, 0.0, 0.0, 0.0),
        NqNumber(acosh(1.0), 0.0, 0.0, 0.0)
    ),
    AcoshTest(
        NqNumber(0.0, 1.0, 0.0, 0.0),
        NqNumber(
            c128Acosh(Complex128(0.0, 1.0)).real,
            c128Acosh(Complex128(0.0, 1.0)).imag,
            0.0,
            0.0
        )
    ),
    AcoshTest(
        NqNumber(0.0, 0.0, 1.0, 0.0),
        NqNumber(
            c128Acosh(Complex128(0.0, 1.0)).real,
            0.0,
            c128Acosh(Complex128(0.0, 1.0)).imag,
            0.0
        )
    ),
    AcoshTest(
        NqNumber(0.0, 0.0, 0.0, 1.0),
        NqNumber(
            c128Acosh(Complex128(0.0, 1.0)).real,
            0.0,
            0.0,
            c128Acosh(Complex128(0.0, 1.0)).imag
        )
    ),
    AcoshTest(
        nqCosh(NqNumber(1.0, 1.0, 1.0, 1.0)),
        NqNumber(1.0, 1.0, 1.0, 1.0)
    ),
    AcoshTest(
        NqNumber(1.0, 1.0, 1.0, 1.0),
        nqLog(
            nqAdd(NqNumber(1.0, 1.0, 1.0, 1.0),
            nqSqrt(nqSub(nqMul(NqNumber(1.0, 1.0, 1.0, 1.0), NqNumber(1.0, 1.0, 1.0, 1.0)), NqNumber(1.0, 0.0, 0.0, 0.0)))))
    )
    // The following fails by a factor of -1.
    // {q: Cosh(NqNumber{1, 1, 1, 1}), want: func() NqNumber {
    // 	q := Cosh(NqNumber{1, 1, 1, 1})
    // 	return Log(Add(q, Sqrt(Sub(Mul(q, q), NqNumber{Real: 1}))))
    // }()},
]

@Test
public func testAcosh(): Unit {
    let tol = 1e-14
    for (test in acoshTests) {
        let got = nqAcosh(test.q)
        @Expect(nqEqualApprox(got, test.want, tol))
    }
}

private class AtanTest {
    let q: NqNumber
    let want: NqNumber
    public init(q: NqNumber, want: NqNumber) {
        this.q = q
        this.want = want
    }
}

private let atanTests: Array<AtanTest> = [
    AtanTest(NqNumber(0.0, 0.0, 0.0, 0.0), NqNumber(0.0, 0.0, 0.0, 0.0)),
    AtanTest(
        NqNumber(1.0, 0.0, 0.0, 0.0),
        NqNumber(Float64.getPI() / 4.0, 0.0, 0.0, 0.0)
    ),
    AtanTest(
        NqNumber(0.0, 0.5, 0.0, 0.0),
        NqNumber(
            c128Atan(Complex128(0.0, 0.5)).real,
            c128Atan(Complex128(0.0, 0.5)).imag,
            0.0,
            0.0
        )
    ),
    AtanTest(
        NqNumber(0.0, 0.0, 0.5, 0.0),
        NqNumber(
            c128Atan(Complex128(0.0, 0.5)).real,
            0.0,
            c128Atan(Complex128(0.0, 0.5)).imag,
            0.0
        )
    ),
    AtanTest(
        NqNumber(0.0, 0.0, 0.0, 0.5),
        NqNumber(
            c128Atan(Complex128(0.0, 0.5)).real,
            0.0,
            0.0,
            c128Atan(Complex128(0.0, 0.5)).imag
        )
    ),
    AtanTest(
        nqTan(NqNumber(1.0, 1.0, 1.0, 1.0)),
        NqNumber(1.0, 1.0, 1.0, 1.0)
    )
]

@Test
public func testAtan(): Unit {
    let tol = 1e-14
    for (test in atanTests) {
        let got = nqAtan(test.q)
        @Expect(nqEqualApprox(got, test.want, tol))
    }
}

private class AtanhTest {
    let q: NqNumber
    let want: NqNumber
    public init(q: NqNumber, want: NqNumber) {
        this.q = q
        this.want = want
    }
}

private let atanhTests: Array<AtanhTest> = [
    AtanhTest(NqNumber(0.0, 0.0, 0.0, 0.0), NqNumber(0.0, 0.0, 0.0, 0.0)),
    AtanhTest(
        NqNumber(1.0, 0.0, 0.0, 0.0),
        NqNumber(nanAtanh(1.0), 0.0, 0.0, 0.0)
    ),
    AtanhTest(
        NqNumber(0.0, 0.5, 0.0, 0.0),
        NqNumber(
            c128Atanh(Complex128(0.0, 0.5)).real,
            c128Atanh(Complex128(0.0, 0.5)).imag,
            0.0,
            0.0
        )
    ),
    AtanhTest(
        NqNumber(0.0, 0.0, 0.5, 0.0),
        NqNumber(
            c128Atanh(Complex128(0.0, 0.5)).real,
            0.0,
            c128Atanh(Complex128(0.0, 0.5)).imag,
            0.0
        )
    ),
    AtanhTest(
        NqNumber(0.0, 0.0, 0.0, 0.5),
        NqNumber(
            c128Atanh(Complex128(0.0, 0.5)).real,
            0.0,
            0.0,
            c128Atanh(Complex128(0.0, 0.5)).imag
        )
    ),
    AtanhTest(
        NqNumber(1.0, 1.0, 1.0, 1.0),
        nqScale(
            0.5,
            nqSub(nqLog(nqAdd(NqNumber(1.0, 0.0, 0.0, 0.0), NqNumber(1.0, 1.0, 1.0, 1.0))),
                nqLog(nqSub(NqNumber(1.0, 0.0, 0.0, 0.0), NqNumber(1.0, 1.0, 1.0, 1.0))))
        )
    ),
    AtanhTest(
        nqTanh(NqNumber(1.0, 0.0, 0.0, 0.0)),
        NqNumber(1.0, 0.0, 0.0, 0.0)
    ),
    AtanhTest(
        nqTanh(NqNumber(0.0, 1.0, 0.0, 0.0)),
        NqNumber(0.0, 1.0, 0.0, 0.0)
    ),
    AtanhTest(
        nqTanh(NqNumber(0.0, 1.0, 1.0, 0.0)),
        NqNumber(0.0, 1.0, 1.0, 0.0)
    ),
    AtanhTest(
        nqTanh(NqNumber(1.0, 1.0, 1.0, 0.0)),
        NqNumber(1.0, 1.0, 1.0, 0.0)
    ),
    // The following fails
    // {q: Tanh(NqNumber{1, 1, 1, 1}), want: NqNumber{1, 1, 1, 1}},
    // but...
    AtanhTest(
        nqTanh(NqNumber(1.0, 1.0, 1.0, 0.0)),
        nqScale(
            0.5,
            nqSub(nqLog(nqAdd(NqNumber(1.0, 0.0, 0.0, 0.0), nqTanh(NqNumber(1.0, 1.0, 1.0, 0.0)))),
                nqLog(nqSub(NqNumber(1.0, 0.0, 0.0, 0.0), nqTanh(NqNumber(1.0, 1.0, 1.0, 0.0)))))
        )
    )
]

@Test
public func testAtanh(): Unit {
    let tol = 1e-14
    for (test in atanhTests) {
        let got = nqAtanh(test.q)
        @Expect(nqEqualApprox(got, test.want, tol))
    }
}
