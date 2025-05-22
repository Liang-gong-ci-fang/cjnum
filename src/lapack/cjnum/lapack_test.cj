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
public func testDlartff(): Unit {
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
@Parallel
public func testIladlr(): Unit {
    @Expect(iladlrTest(lFloat64Implementation))
}

@Test
@Parallel
public func testDlaset(): Unit {
    @Expect(dlasetTest(lFloat64Implementation))
}

@Test
@Parallel
public func testDpttrf(): Unit {
    @Expect(dpttrfTest(lFloat64Implementation))
}

@Test
@Parallel
public func testDlarfg(): Unit {
    @Expect(dlarfgTest(lFloat64Implementation))
}

@Test
@Parallel
public func testDlarf(): Unit {
    @Expect(dlarfTest(lFloat64Implementation))
}

@Test
@Parallel
public func testDlanhs(): Unit {
    @Expect(dlanhsTest(lFloat64Implementation))
}

@Test
@Parallel
public func testDgeql2(): Unit {
    @Expect(dgeql2Test(lFloat64Implementation))
}

@Test
@Parallel
public func testDlanv2(): Unit {
    @Expect(dlanv2Test(lFloat64Implementation))
}

@Test
@Parallel
public func testDgebd2(): Unit {
    @Expect(dgebd2Test(lFloat64Implementation))
}

@Test
@Parallel
public func testDlarfx(): Unit {
    @Expect(dlarfxTest(lFloat64Implementation))
}

@Test
@Parallel
public func testDtrtrs(): Unit {
    @Expect(dtrtrsTest(lFloat64Implementation))
}

@Test
@Skip // 精度
public func testDpttrs(): Unit {
    @Expect(dpttrsTest(lFloat64Implementation))
}

@Test
@Skip // 精度
public func testDlaqr1(): Unit {
    @Expect(dlaqr1Test(lFloat64Implementation))
}

@Test
@Skip // 随机数
public func testDtrtri(): Unit {
    @Expect(dtrtriTest(lFloat64Implementation))
}

@Test
@Skip // 随机数
public func testDgetrf(): Unit {
    @Expect(dgetrfTest(lFloat64Implementation))
}

@Test
@Skip // 随机数
public func testDgecon(): Unit {
    @Expect(dgeconTest(lFloat64Implementation))
}

@Test
@Skip // 随机数
public func testDgetri(): Unit {
    @Expect(dgetriTest(lFloat64Implementation))
}

@Test
@Skip // 随机数
public func testDgetrs(): Unit {
    @Expect(dgetrsTest(lFloat64Implementation))
}
