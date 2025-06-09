package cjnum.num.quat

import cjnum.util.*
import cjnum.complex.*
import std.math.*
import std.unittest.*
import std.unittest.testmacro.*

private class SinTest {
    let q: Number
    let want: Number
    public init(q: Number, want: Number) {
        this.q = q
        this.want = want
    }
}

private let sinTests: Array<SinTest> = [
    SinTest(Number(0.0, 0.0, 0.0, 0.0), Number(0.0, 0.0, 0.0, 0.0)),
    SinTest(Number(Float64.PI / 2.0, 0.0, 0.0, 0.0), Number(1.0, 0.0, 0.0, 0.0)),
    SinTest(
        Number(0.0, Float64.PI / 2.0, 0.0, 0.0),
        Number(0.0, c128Sin(Complex128(0.0, Float64.PI / 2.0)).imag, 0.0, 0.0)
    ),
    SinTest(
        Number(0.0, 0.0, Float64.PI / 2.0, 0.0),
        Number(0.0, 0.0, c128Sin(Complex128(0.0, Float64.PI / 2.0)).imag, 0.0)
    ),
    SinTest(
        Number(0.0, 0.0, 0.0, Float64.PI / 2.0),
        Number(0.0, 0.0, 0.0, c128Sin(Complex128(0.0, Float64.PI / 2.0)).imag)
    ),

    // Exercises from Real Quaternionic Calculus Handbook doi:10.1007/978-3-0348-0622-0
    // Ex 6.159 (a) and (b).
    SinTest(
        Number(1.0, 1.0, 1.0, 1.0),
        Number(sin(1.0) * cosh(sqrt(3.0)), cos(1.0) * sinh(sqrt(3.0)) / sqrt(3.0), cos(1.0) * sinh(sqrt(3.0)) / sqrt(
            3.0), cos(1.0) * sinh(sqrt(3.0)) / sqrt(3.0))
    ),
    SinTest(
        Number(0.0, -2.0, 1.0, 0.0),
        Number(0.0, -2.0 * sinh(sqrt(5.0)) / sqrt(5.0), sinh(sqrt(5.0)) / sqrt(5.0), 0.0)
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
    let q: Number
    let want: Number
    public init(q: Number, want: Number) {
        this.q = q
        this.want = want
    }
}

private let sinhTests: Array<SinhTest> = [
    SinhTest(Number(0.0, 0.0, 0.0, 0.0), Number(0.0, 0.0, 0.0, 0.0)),
    SinhTest(Number(Float64.PI / 2.0, 0.0, 0.0, 0.0), Number(sinh(Float64.PI / 2.0), 0.0, 0.0, 0.0)),
    SinhTest(
        Number(0.0, Float64.PI / 2.0, 0.0, 0.0),
        Number(0.0, c128Sinh(Complex128(0.0, Float64.PI / 2.0)).imag, 0.0, 0.0)
    ),
    SinhTest(
        Number(0.0, 0.0, Float64.PI / 2.0, 0.0),
        Number(0.0, 0.0, c128Sinh(Complex128(0.0, Float64.PI / 2.0)).imag, 0.0)
    ),
    SinhTest(
        Number(0.0, 0.0, 0.0, Float64.PI / 2.0),
        Number(0.0, 0.0, 0.0, c128Sinh(Complex128(0.0, Float64.PI / 2.0)).imag)
    ),
    SinhTest(
        Number(1.0, -1.0, -1.0, 0.0),
        nqScale(
            0.5,
            nqSub(
                nqExp(Number(1.0, -1.0, -1.0, 0.0)),
                nqExp(nqScale(-1.0, Number(1.0, -1.0, -1.0, 0.0)))
            )
        )
    ),
    SinhTest(
        Number(1.0, 1.0, 1.0, 1.0),
        nqScale(
            0.5,
            nqSub(
                nqExp(Number(1.0, 1.0, 1.0, 1.0)),
                nqExp(nqScale(-1.0, Number(1.0, 1.0, 1.0, 1.0)))
            )
        )
    ),
    SinhTest(
        nqAsinh(Number(1.0, 1.0, 1.0, 1.0)),
        Number(1.0, 1.0, 1.0, 1.0)
    ),
    SinhTest(
        nqAsinh(Number(1.0, 1.0, 1.0, 1.0)),
        nqScale(
            0.5,
            nqSub(
                nqExp(nqAsinh(Number(1.0, 1.0, 1.0, 1.0))),
                nqExp(nqScale(-1.0, nqAsinh(Number(1.0, 1.0, 1.0, 1.0))))
            )
        )
    ),
    SinhTest(
        Number(Float64.Inf, 0.0, 0.0, 0.0),
        Number(Float64.Inf, 0.0, 0.0, 0.0)
    ),
    SinhTest(
        Number(Float64.Inf, Float64.PI / 2.0, 0.0, 0.0),
        Number(Float64.Inf, Float64.Inf, 0.0, 0.0)
    ),
    SinhTest(
        Number(Float64.Inf, Float64.PI, 0.0, 0.0),
        Number(-Float64.Inf, Float64.Inf, 0.0, 0.0)
    ),
    SinhTest(
        Number(Float64.Inf, 3.0 * Float64.PI / 2.0, 0.0, 0.0),
        Number(-Float64.Inf, -Float64.Inf, 0.0, 0.0)
    ),
    SinhTest(
        Number(Float64.Inf, 2.0 * Float64.PI, 0.0, 0.0),
        Number(Float64.Inf, -Float64.Inf, 0.0, 0.0)
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
    let q: Number
    let want: Number
    public init(q: Number, want: Number) {
        this.q = q
        this.want = want
    }
}

private let cosTests: Array<CosTest> = [
    CosTest(Number(0.0, 0.0, 0.0, 0.0), Number(1.0, 0.0, 0.0, 0.0)),
    CosTest(
        Number(Float64.PI / 2.0, 0.0, 0.0, 0.0),
        Number(0.0, 0.0, 0.0, 0.0)
    ),
    CosTest(
        Number(0.0, Float64.PI / 2.0, 0.0, 0.0),
        Number(c128Cos(Complex128(0.0, Float64.PI / 2.0)).real, 0.0, 0.0, 0.0)
    ),
    CosTest(
        Number(0.0, 0.0, Float64.PI / 2.0, 0.0),
        Number((c128Cos(Complex128(0.0, Float64.PI / 2.0))).real, 0.0, 0.0, 0.0)
    ),
    CosTest(
        Number(0.0, 0.0, 0.0, Float64.PI / 2.0),
        Number((c128Cos(Complex128(0.0, Float64.PI / 2.0))).real, 0.0, 0.0, 0.0)
    ),
    CosTest(
        Number(1.0, 1.0, 1.0, 1.0),
        Number(
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
    let q: Number
    let want: Number
    public init(q: Number, want: Number) {
        this.q = q
        this.want = want
    }
}

private let coshTests: Array<CoshTest> = [
    CoshTest(Number(0.0, 0.0, 0.0, 0.0), Number(1.0, 0.0, 0.0, 0.0)),
    CoshTest(
        Number(Float64.PI / 2.0, 0.0, 0.0, 0.0),
        Number(cosh(Float64.PI / 2.0), 0.0, 0.0, 0.0)
    ),
    CoshTest(
        Number(0.0, Float64.PI / 2.0, 0.0, 0.0),
        Number(0.0, c128Cosh(Complex128(0.0, Float64.PI / 2.0)).imag, 0.0, 0.0)
    ),
    CoshTest(
        Number(0.0, 0.0, Float64.PI / 2.0, 0.0),
        Number(0.0, 0.0, c128Cosh(Complex128(0.0, Float64.PI / 2.0)).imag, 0.0)
    ),
    CoshTest(
        Number(0.0, 0.0, 0.0, Float64.PI / 2.0),
        Number(0.0, 0.0, 0.0, c128Cosh(Complex128(0.0, Float64.PI / 2.0)).imag)
    ),
    CoshTest(
        Number(1.0, -1.0, -1.0, 0.0),
        nqScale(0.5, nqAdd(nqExp(Number(1.0, -1.0, -1.0, 0.0)), nqExp(nqScale(-1.0, Number(1.0, -1.0, -1.0, 0.0)))))
    ),
    CoshTest(
        Number(1.0, 1.0, 1.0, 1.0),
        nqScale(0.5, nqAdd(nqExp(Number(1.0, 1.0, 1.0, 1.0)), nqExp(nqScale(-1.0, Number(1.0, 1.0, 1.0, 1.0)))))
    ),
    CoshTest(
        Number(Float64.Inf, 0.0, 0.0, 0.0),
        Number(Float64.Inf, 0.0, 0.0, 0.0)
    ),
    CoshTest(
        Number(Float64.Inf, Float64.PI / 2.0, 0.0, 0.0),
        Number(Float64.Inf, Float64.Inf, 0.0, 0.0)
    ),
    CoshTest(
        Number(Float64.Inf, Float64.PI, 0.0, 0.0),
        Number(-Float64.Inf, Float64.Inf, 0.0, 0.0)
    ),
    CoshTest(
        Number(Float64.Inf, 3.0 * Float64.PI / 2.0, 0.0, 0.0),
        Number(-Float64.Inf, -Float64.Inf, 0.0, 0.0)
    ),
    CoshTest(
        Number(Float64.Inf, 2.0 * Float64.PI, 0.0, 0.0),
        Number(Float64.Inf, -Float64.Inf, 0.0, 0.0)
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
    let q: Number
    let want: Number
    public init(q: Number, want: Number) {
        this.q = q
        this.want = want
    }
}

private let tanTests: Array<TanTest> = [
    TanTest(Number(0.0, 0.0, 0.0, 0.0), Number(1.0, 0.0, 0.0, 0.0)),
    TanTest(
        Number(Float64.PI / 4.0, 0.0, 0.0, 0.0),
        Number(tan(Float64.PI / 4.0), 0.0, 0.0, 0.0)
    ),
    TanTest(
        Number(0.0, Float64.PI / 4.0, 0.0, 0.0),
        Number(0.0, c128Tan(Complex128(0.0, Float64.PI / 4.0)).imag, 0.0, 0.0)
    ),
    TanTest(
        Number(0.0, 0.0, Float64.PI / 4.0, 0.0),
        Number(0.0, 0.0, c128Tan(Complex128(0.0, Float64.PI / 4.0)).imag, 0.0)
    ),
    TanTest(
        Number(0.0, 0.0, 0.0, Float64.PI / 4.0),
        Number(0.0, 0.0, 0.0, c128Tan(Complex128(0.0, Float64.PI / 4.0)).imag)
    ),
    TanTest(
        Number(0.0, 1.0, 0.0, 0.0),
        nqMul(
            nqSin(Number(0.0, 1.0, 0.0, 0.0)),
            nqInv(nqCos(Number(0.0, 1.0, 0.0, 0.0)))
        )
    ),
    TanTest(
        Number(1.0, 1.0, 1.0, 1.0),
        nqMul(
            nqSin(Number(1.0, 1.0, 1.0, 1.0)),
            nqInv(nqCos(Number(1.0, 1.0, 1.0, 1.0)))
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
    let q: Number
    let want: Number
    public init(q: Number, want: Number) {
        this.q = q
        this.want = want
    }
}

private let tanhTests: Array<TanhTest> = [
    TanhTest(Number(0.0, 0.0, 0.0, 0.0), Number(0.0, 0.0, 0.0, 0.0)),
    TanhTest(
        Number(Float64.PI / 4.0, 0.0, 0.0, 0.0),
        Number(tanh(Float64.PI / 4.0), 0.0, 0.0, 0.0)
    ),
    TanhTest(
        Number(0.0, Float64.PI / 4.0, 0.0, 0.0),
        Number(0.0, c128Tanh(Complex128(0.0, Float64.PI / 4.0)).imag, 0.0, 0.0)
    ),
    TanhTest(
        Number(0.0, 0.0, Float64.PI / 4.0, 0.0),
        Number(0.0, 0.0, c128Tanh(Complex128(0.0, Float64.PI / 4.0)).imag, 0.0)
    ),
    TanhTest(
        Number(0.0, 0.0, 0.0, Float64.PI / 4.0),
        Number(0.0, 0.0, 0.0, c128Tanh(Complex128(0.0, Float64.PI / 4.0)).imag)
    ),
    TanhTest(
        Number(0.0, 1.0, 0.0, 0.0),
        nqMul(nqSinh(Number(0.0, 1.0, 0.0, 0.0)), nqInv(nqCosh(Number(0.0, 1.0, 0.0, 0.0))))
    ),
    TanhTest(
        Number(1.0, 1.0, 1.0, 1.0),
        nqMul(nqSinh(Number(1.0, 1.0, 1.0, 1.0)), nqInv(nqCosh(Number(1.0, 1.0, 1.0, 1.0))))
    ),
    TanhTest(
        Number(Float64.Inf, 0.0, 0.0, 0.0),
        Number(1.0, 0.0, 0.0, 0.0)
    ),
    TanhTest(
        Number(Float64.Inf, Float64.PI / 4.0, 0.0, 0.0),
        Number(1.0, 0.0 * sin(Float64.PI / 2.0), 0.0, 0.0)
    ),
    TanhTest(
        Number(Float64.Inf, Float64.PI / 2.0, 0.0, 0.0),
        Number(1.0, 0.0 * sin(Float64.PI), 0.0, 0.0)
    ),
    TanhTest(
        Number(Float64.Inf, 3.0 * Float64.PI / 4.0, 0.0, 0.0),
        Number(1.0, 0.0 * sin(3.0 * Float64.PI / 2.0), 0.0, 0.0)
    ),
    TanhTest(
        Number(Float64.Inf, Float64.PI, 0.0, 0.0),
        Number(1.0, 0.0 * sin(2.0 * Float64.PI), 0.0, 0.0)
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
    let q: Number
    let want: Number
    public init(q: Number, want: Number) {
        this.q = q
        this.want = want
    }
}

private let asinTests: Array<AsinTest> = [
    AsinTest(Number(0.0, 0.0, 0.0, 0.0), Number(0.0, 0.0, 0.0, 0.0)),
    AsinTest(
        Number(1.0, 0.0, 0.0, 0.0),
        Number(Float64.PI / 2.0, 0.0, 0.0, 0.0)
    ),
    AsinTest(
        Number(0.0, 1.0, 0.0, 0.0),
        Number(0.0, c128Asinh(Complex128(1.0, 0.0)).real, 0.0, 0.0)
    ),
    AsinTest(
        Number(0.0, 0.0, 1.0, 0.0),
        Number(0.0, 0.0, c128Asinh(Complex128(1.0, 0.0)).real, 0.0)
    ),
    AsinTest(
        Number(0.0, 0.0, 0.0, 1.0),
        Number(0.0, 0.0, 0.0, c128Asinh(Complex128(1.0, 0.0)).real)
    ),
    AsinTest(
        nqSin(Number(1.0, 1.0, 1.0, 1.0)),
        Number(1.0, 1.0, 1.0, 1.0)
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
    let q: Number
    let want: Number
    public init(q: Number, want: Number) {
        this.q = q
        this.want = want
    }
}

private let asinhTests: Array<AsinhTest> = [
    AsinhTest(Number(0.0, 0.0, 0.0, 0.0), Number(0.0, 0.0, 0.0, 0.0)),
    AsinhTest(
        Number(1.0, 0.0, 0.0, 0.0),
        Number(asinh(1.0), 0.0, 0.0, 0.0)
    ),
    AsinhTest(
        Number(0.0, 1.0, 0.0, 0.0),
        Number(0.0, Float64.PI / 2.0, 0.0, 0.0)
    ),
    AsinhTest(
        Number(0.0, 0.0, 1.0, 0.0),
        Number(0.0, 0.0, Float64.PI / 2.0, 0.0)
    ),
    AsinhTest(
        Number(0.0, 0.0, 0.0, 1.0),
        Number(0.0, 0.0, 0.0, Float64.PI / 2.0)
    ),
    AsinhTest(
        Number(1.0, 1.0, 1.0, 1.0),
        nqLog(
            nqAdd(
                Number(1.0, 1.0, 1.0, 1.0),
                nqSqrt(
                    nqAdd(
                        nqMul(Number(1.0, 1.0, 1.0, 1.0), Number(1.0, 1.0, 1.0, 1.0)),
                        Number(1.0, 0.0, 0.0, 0.0)
                    )
                )
            )
        )
    ),
    AsinhTest(
        nqSinh(Number(1.0, 0.0, 0.0, 0.0)),
        Number(1.0, 0.0, 0.0, 0.0)
    ),
    AsinhTest(
        nqSinh(Number(0.0, 1.0, 0.0, 0.0)),
        Number(0.0, 1.0, 0.0, 0.0)
    ),
    AsinhTest(
        nqSinh(Number(0.0, 1.0, 1.0, 0.0)),
        Number(0.0, 1.0, 1.0, 0.0)
    ),
    AsinhTest(
        nqSinh(Number(1.0, 1.0, 1.0, 0.0)),
        Number(1.0, 1.0, 1.0, 0.0)
    ),
    // The following fails:
    // {q: Sinh(Number{1, 1, 1, 1}), want: Number{1, 1, 1, 1}},
    // but this passes...
    AsinhTest(
        nqSinh(Number(1.0, 1.0, 1.0, 1.0)),
        nqLog(
            nqAdd(
                nqSinh(Number(1.0, 1.0, 1.0, 1.0)),
                nqSqrt(
                    nqAdd(
                        nqMul(nqSinh(Number(1.0, 1.0, 1.0, 1.0)), nqSinh(Number(1.0, 1.0, 1.0, 1.0))),
                        Number(1.0, 0.0, 0.0, 0.0)
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
    let q: Number
    let want: Number
    public init(q: Number, want: Number) {
        this.q = q
        this.want = want
    }
}

private let acosTests: Array<AcosTest> = [
    AcosTest(Number(0.0, 0.0, 0.0, 0.0), Number(Float64.PI / 2.0, 0.0, 0.0, 0.0)),
    AcosTest(
        Number(1.0, 0.0, 0.0, 0.0),
        Number(0.0, 0.0, 0.0, 0.0)
    ),
    AcosTest(
        Number(0.0, 1.0, 0.0, 0.0),
        Number(
            c128Acos(Complex128(0.0, 1.0)).real,
            c128Acos(Complex128(0.0, 1.0)).imag,
            0.0,
            0.0
        )
    ),
    AcosTest(
        Number(0.0, 0.0, 1.0, 0.0),
        Number(
            c128Acos(Complex128(0.0, 1.0)).real,
            0.0,
            c128Acos(Complex128(0.0, 1.0)).imag,
            0.0
        )
    ),
    AcosTest(
        Number(0.0, 0.0, 0.0, 1.0),
        Number(
            c128Acos(Complex128(0.0, 1.0)).real,
            0.0,
            0.0,
            c128Acos(Complex128(0.0, 1.0)).imag
        )
    ),
    AcosTest(
        nqCos(Number(1.0, 1.0, 1.0, 1.0)),
        Number(1.0, 1.0, 1.0, 1.0)
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
    let q: Number
    let want: Number
    public init(q: Number, want: Number) {
        this.q = q
        this.want = want
    }
}

private let acoshTests: Array<AcoshTest> = [
    AcoshTest(Number(0.0, 0.0, 0.0, 0.0), Number(Float64.PI / 2.0, 0.0, 0.0, 0.0)),
    AcoshTest(
        Number(1.0, 0.0, 0.0, 0.0),
        Number(acosh(1.0), 0.0, 0.0, 0.0)
    ),
    AcoshTest(
        Number(0.0, 1.0, 0.0, 0.0),
        Number(
            c128Acosh(Complex128(0.0, 1.0)).real,
            c128Acosh(Complex128(0.0, 1.0)).imag,
            0.0,
            0.0
        )
    ),
    AcoshTest(
        Number(0.0, 0.0, 1.0, 0.0),
        Number(
            c128Acosh(Complex128(0.0, 1.0)).real,
            0.0,
            c128Acosh(Complex128(0.0, 1.0)).imag,
            0.0
        )
    ),
    AcoshTest(
        Number(0.0, 0.0, 0.0, 1.0),
        Number(
            c128Acosh(Complex128(0.0, 1.0)).real,
            0.0,
            0.0,
            c128Acosh(Complex128(0.0, 1.0)).imag
        )
    ),
    AcoshTest(
        nqCosh(Number(1.0, 1.0, 1.0, 1.0)),
        Number(1.0, 1.0, 1.0, 1.0)
    ),
    AcoshTest(
        Number(1.0, 1.0, 1.0, 1.0),
        nqLog(
            nqAdd(Number(1.0, 1.0, 1.0, 1.0),
            nqSqrt(nqSub(nqMul(Number(1.0, 1.0, 1.0, 1.0), Number(1.0, 1.0, 1.0, 1.0)), Number(1.0, 0.0, 0.0, 0.0)))))
    )
    // The following fails by a factor of -1.
    // {q: Cosh(Number{1, 1, 1, 1}), want: func() Number {
    // 	q := Cosh(Number{1, 1, 1, 1})
    // 	return Log(Add(q, Sqrt(Sub(Mul(q, q), Number{Real: 1}))))
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
    let q: Number
    let want: Number
    public init(q: Number, want: Number) {
        this.q = q
        this.want = want
    }
}

private let atanTests: Array<AtanTest> = [
    AtanTest(Number(0.0, 0.0, 0.0, 0.0), Number(0.0, 0.0, 0.0, 0.0)),
    AtanTest(
        Number(1.0, 0.0, 0.0, 0.0),
        Number(Float64.PI / 4.0, 0.0, 0.0, 0.0)
    ),
    AtanTest(
        Number(0.0, 0.5, 0.0, 0.0),
        Number(
            c128Atan(Complex128(0.0, 0.5)).real,
            c128Atan(Complex128(0.0, 0.5)).imag,
            0.0,
            0.0
        )
    ),
    AtanTest(
        Number(0.0, 0.0, 0.5, 0.0),
        Number(
            c128Atan(Complex128(0.0, 0.5)).real,
            0.0,
            c128Atan(Complex128(0.0, 0.5)).imag,
            0.0
        )
    ),
    AtanTest(
        Number(0.0, 0.0, 0.0, 0.5),
        Number(
            c128Atan(Complex128(0.0, 0.5)).real,
            0.0,
            0.0,
            c128Atan(Complex128(0.0, 0.5)).imag
        )
    ),
    AtanTest(
        nqTan(Number(1.0, 1.0, 1.0, 1.0)),
        Number(1.0, 1.0, 1.0, 1.0)
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
    let q: Number
    let want: Number
    public init(q: Number, want: Number) {
        this.q = q
        this.want = want
    }
}

private let atanhTests: Array<AtanhTest> = [
    AtanhTest(Number(0.0, 0.0, 0.0, 0.0), Number(0.0, 0.0, 0.0, 0.0)),
    AtanhTest(
        Number(1.0, 0.0, 0.0, 0.0),
        Number(nanAtanh(1.0), 0.0, 0.0, 0.0)
    ),
    AtanhTest(
        Number(0.0, 0.5, 0.0, 0.0),
        Number(
            c128Atanh(Complex128(0.0, 0.5)).real,
            c128Atanh(Complex128(0.0, 0.5)).imag,
            0.0,
            0.0
        )
    ),
    AtanhTest(
        Number(0.0, 0.0, 0.5, 0.0),
        Number(
            c128Atanh(Complex128(0.0, 0.5)).real,
            0.0,
            c128Atanh(Complex128(0.0, 0.5)).imag,
            0.0
        )
    ),
    AtanhTest(
        Number(0.0, 0.0, 0.0, 0.5),
        Number(
            c128Atanh(Complex128(0.0, 0.5)).real,
            0.0,
            0.0,
            c128Atanh(Complex128(0.0, 0.5)).imag
        )
    ),
    AtanhTest(
        Number(1.0, 1.0, 1.0, 1.0),
        nqScale(
            0.5,
            nqSub(nqLog(nqAdd(Number(1.0, 0.0, 0.0, 0.0), Number(1.0, 1.0, 1.0, 1.0))),
                nqLog(nqSub(Number(1.0, 0.0, 0.0, 0.0), Number(1.0, 1.0, 1.0, 1.0))))
        )
    ),
    AtanhTest(
        nqTanh(Number(1.0, 0.0, 0.0, 0.0)),
        Number(1.0, 0.0, 0.0, 0.0)
    ),
    AtanhTest(
        nqTanh(Number(0.0, 1.0, 0.0, 0.0)),
        Number(0.0, 1.0, 0.0, 0.0)
    ),
    AtanhTest(
        nqTanh(Number(0.0, 1.0, 1.0, 0.0)),
        Number(0.0, 1.0, 1.0, 0.0)
    ),
    AtanhTest(
        nqTanh(Number(1.0, 1.0, 1.0, 0.0)),
        Number(1.0, 1.0, 1.0, 0.0)
    ),
    // The following fails
    // {q: Tanh(Number{1, 1, 1, 1}), want: Number{1, 1, 1, 1}},
    // but...
    AtanhTest(
        nqTanh(Number(1.0, 1.0, 1.0, 0.0)),
        nqScale(
            0.5,
            nqSub(nqLog(nqAdd(Number(1.0, 0.0, 0.0, 0.0), nqTanh(Number(1.0, 1.0, 1.0, 0.0)))),
                nqLog(nqSub(Number(1.0, 0.0, 0.0, 0.0), nqTanh(Number(1.0, 1.0, 1.0, 0.0)))))
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
