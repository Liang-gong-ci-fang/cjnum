package cjnum.lapack.cjnum

import std.unittest.*
import std.unittest.testmacro.*
import cjnum.lapack.testlapack.*

@Test
@Parallel
public func testDrscl(): Unit {
    @Expect(drsclTest(lFloat64Implementation))
}

@Test
@Parallel
public func testDlassq(): Unit {
    @Expect(dlassqTest(lFloat64Implementation))
}

@Test
@Parallel
public func testDlange(): Unit {
    @Expect(dlangeTest(lFloat64Implementation))
}

@Test
@Parallel
public func testDlatrs(): Unit {
    @Expect(dlatrsTest(lFloat64Implementation))
}

@Test
@Parallel
public func testDlacn2(): Unit {
    @Expect(dlacn2Test(lFloat64Implementation))
}

@Test
@Parallel
public func testDtrti2(): Unit {
    @Expect(dtrti2Test(lFloat64Implementation))
}

@Test
@Parallel
public func testDgetf2(): Unit {
    @Expect(dgetf2Test(lFloat64Implementation))
}

@Test
@Parallel
public func testDlaswp(): Unit {
    @Expect(dlaswpTest(lFloat64Implementation))
}

@Test
@Parallel
public func testDlartf(): Unit {
    @Expect(dlartgTest(lFloat64Implementation))
}

@Test
@Parallel
public func testDlas2(): Unit {
    @Expect(dlas2Test(lFloat64Implementation))
}

@Test
@Parallel
public func testDlae2(): Unit {
    @Expect(dlae2Test(lFloat64Implementation))
}

@Test
@Parallel
public func testDlacpy(): Unit {
    @Expect(dlacpyTest(lFloat64Implementation))
}

@Test
@Parallel
public func testDlansy(): Unit {
    @Expect(dlansyTest(lFloat64Implementation))
}

@Test
@Parallel
public func testIladlc(): Unit {
    @Expect(iladlcTest(lFloat64Implementation))
}

@Test
@Skip
public func testDtrtri(): Unit {
    @Expect(dtrtriTest(lFloat64Implementation))
}

@Test
@Skip
public func testDgetrf(): Unit {
    @Expect(dgetrfTest(lFloat64Implementation))
}

@Test
@Skip
public func testDgecon(): Unit {
    @Expect(dgeconTest(lFloat64Implementation))
}

@Test
@Skip
public func testDgetri(): Unit {
    @Expect(dgetriTest(lFloat64Implementation))
}

@Test
@Skip
public func testDgetrs(): Unit {
    @Expect(dgetrsTest(lFloat64Implementation))
}