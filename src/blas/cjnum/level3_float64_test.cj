package cjnum.blas.cjnum

import cjnum.blas.testblas.*
import std.unittest.*
import std.unittest.testmacro.*

@Test
public func testDtrsm(): Unit {
    @Expect(dtrsmTest(nFloat64Implementation))
}

@Test
public func testDsymm(): Unit {
    @Expect(dsymmTest(nFloat64Implementation))
}

@Test
public func testDsyrk(): Unit {
    @Expect(dsyrkTest(nFloat64Implementation))
}

@Test
public func testDsyr2k(): Unit {
    @Expect(dsyr2kTest(nFloat64Implementation))
}

@Test
public func testDtrmm(): Unit {
    @Expect(dtrmmTest(nFloat64Implementation))
}

@Test
public func testDgemm(): Unit {
    @Expect(dgemmTest(nFloat64Implementation))
}