package cjnum.num.quat

import cjnum.floats.scalar.*
import std.math.*
import std.unittest.*
import std.unittest.testmacro.*

private class InvTest {
    let q: NqNumber
    let wantNaN: Bool
    public init(q: NqNumber, wantNaN: Bool) {
        this.q = q
        this.wantNaN = wantNaN
    }
}

/*
因为仓颉不允许传入空值，所以这里依据原项目要求(对于大多数测试用例，wantNaN 被设置为默认值 false，
表示期望得到正常的四元数结果)在测试用例中以false代替了空值wantNaN
 */
private let invTests: Array<InvTest> = [
    InvTest(NqNumber(1.0, 1.0, 1.0, 1.0), false),
    InvTest(NqNumber(3.0, -1.0, 5.0, -40.0), false),
    InvTest(NqNumber(1e6, -1e5, 4.0, -10.0), false),
    InvTest(NqNumber(0.0, 1.0, 1.0, 1.0), false),
    InvTest(NqNumber(1.0, 0.0, 1.0, 1.0), false),
    InvTest(NqNumber(1.0, 1.0, 0.0, 1.0), false),
    InvTest(NqNumber(1.0, 1.0, 1.0, 0.0), false),
    InvTest(NqNumber(0.0, 0.0, 0.0, 0.0), true)
]

@Test
public func testInv(): Unit {
    let tol = 1e-14
    for (test in invTests) {
        let got = nqMul(test.q, nqInv(test.q))
        if (test.wantNaN) {
            @Expect(nqIsNaN(got))
            continue
        }
        @Expect(equalWithinAbsOrRel(got.real, 1.0, tol, tol) && equalWithinAbsOrRel(nqAbs(got), 1.0, tol, tol))
    }
}
