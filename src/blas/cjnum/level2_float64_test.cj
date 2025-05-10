package cjnum.blas.cjnum

import cjnum.blas.testblas.*
import std.unittest.*
import std.unittest.testmacro.*


@Test
public func testDger(): Unit {
    @Expect(dgerTest(nFloat64Implementation))
}

@Test func testDgbmv(): Unit {
    @Expect(dgbmvTest(nFloat64Implementation))
}

@Test
public func testDgemv(): Unit {
    @Expect(dgemvTest(nFloat64Implementation))
}

@Test
public func testDtrmv(): Unit {
    @Expect(dtrmvTest(nFloat64Implementation))
}

@Test
public func testDtrsv(): Unit {
    @Expect(dtrsvTest(nFloat64Implementation))
}

@Test
public func testDsymv(): Unit {
    @Expect(dsymvTest(nFloat64Implementation))
}

@Test
public func testDtbmv(): Unit {
    @Expect(dtbmvTest(nFloat64Implementation))
}

@Test
public func testDtpmv(): Unit {
    @Expect(dtpmvTest(nFloat64Implementation))
}

@Test
public func testDtbsv(): Unit {
    @Expect(dtbsvTest(nFloat64Implementation))
}

@Test
public func testDsbmv(): Unit {
    @Expect(dsbmvTest(nFloat64Implementation))
}

@Test
public func testDsyr(): Unit {
    @Expect(dsyrTest(nFloat64Implementation))
}

@Test
public func testDsyr2(): Unit {
    @Expect(dsyr2Test(nFloat64Implementation))
}

@Test
public func testDtpsv(): Unit {
    @Expect(dtpsvTest(nFloat64Implementation))
}

@Test
public func testDspmv(): Unit {
    @Expect(dspmvTest(nFloat64Implementation))
}

@Test
public func testDspr(): Unit {
    @Expect(dsprTest(nFloat64Implementation))
}

@Test
public func testDspr2(): Unit {
    @Expect(dspr2Test(nFloat64Implementation))
}

@Test
public func testDtxmv(): Unit {
    @Expect(dtxmvTest(nFloat64Implementation))
}