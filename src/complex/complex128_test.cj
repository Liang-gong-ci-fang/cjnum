package cjnum.complex

import std.unittest.*
import std.unittest.testmacro.*

@Test
public class Complex128Test {
    @TestCase
    public func testIsNaN(): Unit {
        let c1 = Complex128(0.0, 0.0)
        @Expect(!c1.isNaN())
        let c2 = Complex128.NaN
        @Expect(c2.isNaN())
        let c3 = Complex128(0.0, Float64.NaN)
        @Expect(c3.isNaN())
        let c4 = Complex128(Float64.NaN, 0.0)
        @Expect(c4.isNaN())
    }

    @TestCase
    public func testEqualAndNotEqual(): Unit {
        let c1 = Complex128(1.0, 2.0)
        let c2 = Complex128(1.0, 2.0)
        @Expect(c1 == c2)
        @Expect(!(c1!= c2))
        let c3 = Complex128(2.0, 1.0)
        @Expect(!(c1 == c3))
        @Expect(c1 != c3)
    }
}
