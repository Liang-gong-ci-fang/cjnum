package cjnum.blas.cjnum

import std.unittest.*
import std.unittest.testmacro.*
import cjnum.blas.testblas.*

let float64Implementation = Float64Implementation()

@Test
public func testDnrm2() {
    @Expect(dnrm2Test(float64Implementation))
}

@Test
public func testDasum() {
    @Expect(dasumTest(float64Implementation))
}

@Test
public func testIdamax() {
    @Expect(idamaxTest(float64Implementation))
}

@Test
public func testDswap() {
    @Expect(dswapTest(float64Implementation))
}

@Test
public func testDcopy() {
    @Expect(dcopyTest(float64Implementation))
}

@Test
public func testDaxpy() {
    @Expect(daxpyTest(float64Implementation))
}

@Test
public func testDrotg() {
    @Expect(drotgTest(float64Implementation, false))
}

@Test
public func testDrotmg() {
    @Expect(drotmgTest(float64Implementation))
}

@Test
public func testDrot() {
    @Expect(drotTest(float64Implementation))
}

@Test
public func testDrotm() {
    @Expect(drotmTest(float64Implementation))
}

@Test
public func testDscal() {
    @Expect(dscalTest(float64Implementation))
}