package cjnum.blas.cjnum

import std.unittest.*
import std.unittest.testmacro.*
import cjnum.blas.testblas.*

let nFloat64Implementation = NFloat64Implementation()

@Test
public func testDnrm2() {
    @Expect(dnrm2Test(nFloat64Implementation))
}

@Test
public func testDasum() {
    @Expect(dasumTest(nFloat64Implementation))
}

@Test
public func testIdamax() {
    @Expect(idamaxTest(nFloat64Implementation))
}

@Test
public func testDswap() {
    @Expect(dswapTest(nFloat64Implementation))
}

@Test
public func testDcopy() {
    @Expect(dcopyTest(nFloat64Implementation))
}

@Test
public func testDaxpy() {
    @Expect(daxpyTest(nFloat64Implementation))
}

@Test
public func testDrotg() {
    @Expect(drotgTest(nFloat64Implementation, false))
}

@Test
public func testDrotmg() {
    @Expect(drotmgTest(nFloat64Implementation))
}

@Test
public func testDrot() {
    @Expect(drotTest(nFloat64Implementation))
}

@Test
public func testDrotm() {
    @Expect(drotmTest(nFloat64Implementation))
}

@Test
public func testDscal() {
    @Expect(dscalTest(nFloat64Implementation))
}