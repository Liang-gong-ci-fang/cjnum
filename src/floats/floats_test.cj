package cjnum.floats

import std.math.*
import std.unittest.*
import std.unittest.testmacro.*

// TODO

@Test
public func testEquals() {
    var s1 = [1.0, 2.0, 3.0, 4.0]
    var s2 = [1.0, 2.0, 3.0, 4.0]
    @Expect(equal(s1, s2))
    s2 = [1.0, 2.0, 3.0, 4.0 + 1e-14]
    @Expect(!equal(s1, s2))
    @Expect(!equal(s1, []))
}

@Test
public func testSame() {
    var s1 = [1.0, 2.0, 3.0, 4.0]
    var s2 = [1.0, 2.0, 3.0, 4.0]
    @Expect(same(s1, s2))
    s2 = [1.0, 2.0, 3.0, 4.0 + 1e-14]
    @Expect(!same(s1, s2))
    @Expect(!same(s1, []))
    s1 = [1.0, 2.0, Float64.NaN, 4.0]
    s2 = [1.0, 2.0, Float64.NaN, 4.0]
    @Expect(same(s1, s2))
    s1 = [1.0, 2.0, Float64.NaN, 4.0]
    s2 = [1.0, Float64.NaN, 3.0, 4.0]
    @Expect(!same(s1, s2))
} 

// TODO