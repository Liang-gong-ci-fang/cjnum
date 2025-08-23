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
@Parallel
public func testDlansb(): Unit {
    @Expect(dlansbTest(lFloat64Implementation))
}

@Test
@Parallel
public func testDlanst(): Unit {
    @Expect(dlanstTest(lFloat64Implementation))
}

@Test
@Parallel
public func testDlantb(): Unit {
    @Expect(dlantbTest(lFloat64Implementation))
}

@Test
@Parallel
public func testDlantr(): Unit {
    @Expect(dlantrTest(lFloat64Implementation))
}

@Test
@Parallel
public func testDlasr(): Unit {
    @Expect(dlasrTest(lFloat64Implementation))
}

@Test
@Parallel
public func testDlays2(): Unit {
    @Expect(dlasy2Test(lFloat64Implementation))
}

@Test
@Parallel
public func testDlascl(): Unit {
    @Expect(dlasclTest(lFloat64Implementation))
}

@Test
@Parallel
public func testDlaexc(): Unit {
    @Expect(dlaexcTest(lFloat64Implementation))
}

@Test
@Parallel
public func testDlasv2(): Unit {
    @Expect(dlasv2Test(lFloat64Implementation))
}

@Test
@Parallel
public func testDlauu2(): Unit {
    @Expect(dlauu2Test(lFloat64Implementation))
}

@Test
@Parallel
public func testDpotf2(): Unit {
    @Expect(dpotf2Test(lFloat64Implementation))
}

@Test
@Parallel
public func testdpbtf2(): Unit {
    @Expect(dpbtf2Test(lFloat64Implementation))
}

@Test
@Parallel
public func testdgetc2(): Unit {
    @Expect(dgetc2Test(lFloat64Implementation))
}

@Test
@Parallel
public func testdgesc2(): Unit {
    @Expect(dgesc2Test(lFloat64Implementation))
}

@Test
@Parallel
public func testDgghrd(): Unit {
    @Expect(dgghrdTest(lFloat64Implementation))
}

@Test
@Parallel
public func testDgtsv(): Unit {
    @Expect(dgtsvTest(lFloat64Implementation))
}

@Test
@Parallel
public func testDgehd2(): Unit {
    @Expect(dgehd2Test(lFloat64Implementation))
}

@Test
@Parallel
public func testDlatrd(): Unit {
    @Expect(dlatrdTest(lFloat64Implementation))
}

@Test
@Parallel
public func testDlahqr(): Unit {
    @Expect(dlahqrTest(lFloat64Implementation))
}

@Test
@Parallel
public func testDgeqr2(): Unit {
    @Expect(dgeqr2Test(lFloat64Implementation))
}

@Test
@Parallel
public func testDlaqr5(): Unit {
    @Expect(dlaqr5Test(lFloat64Implementation))
}

@Test
@Parallel
public func testDlarft(): Unit {
    @Expect(dlarftTest(lFloat64Implementation))
}

@Test
@Parallel
public func testDorg2l(): Unit {
    @Expect(dorg2lTest(lFloat64Implementation))
}

@Test
@Parallel
public func testDtrexc(): Unit {
    @Expect(dtrexcTest(lFloat64Implementation))
}

@Test
@Parallel
public func testDlasrt(): Unit {
    @Expect(dlasrtTest(lFloat64Implementation))
}

@Test
@Parallel
public func testDlaev2(): Unit {
    @Expect(dlaev2Test(lFloat64Implementation))
}

@Test
@Parallel
public func testDlarfb(): Unit {
    @Expect(dlarfbTest(lFloat64Implementation))
}

@Test
@Parallel
public func testDpbtrf(): Unit {
    @Expect(dpbtrfTest(lFloat64Implementation))
}

@Test
@Parallel
public func testDpstf2(): Unit {
    @Expect(dpstf2Test(lFloat64Implementation))
}

@Test
@Parallel
public func testDlahr2(): Unit {
    @Expect(dlahr2Test(lFloat64Implementation))
}

@Test
@Parallel
public func testDlag2(): Unit {
    @Expect(dlag2Test(lFloat64Implementation))
}

@Test
@Parallel
public func testDlags2(): Unit {
    @Expect(dlags2Test(lFloat64Implementation))
}

@Test
@Parallel
public func testDlagtm(): Unit {
    @Expect(dlagtmTest(lFloat64Implementation))
}

@Test
@Parallel
public func testDlangt(): Unit {
    @Expect(dlangtTest(lFloat64Implementation))
}

@Test
@Parallel
public func testDlaln2(): Unit {
    @Expect(dlaln2Test(lFloat64Implementation))
}

@Test
@Parallel
public func testDgebak(): Unit {
    @Expect(dgebakTest(lFloat64Implementation))
}

@Test
@Parallel
public func testDgelq2(): Unit {
    @Expect(dgelq2Test(lFloat64Implementation))
}

@Test
@Parallel
public func testDgerq2(): Unit {
    @Expect(dgerq2Test(lFloat64Implementation))
}

@Test
@Parallel
public func testDlatbs(): Unit {
    @Expect(dlatbsTest(lFloat64Implementation))
}

@Test
@Parallel
public func testDlapmr(): Unit {
    @Expect(dlapmrTest(lFloat64Implementation))
}

@Test
@Parallel
public func testDlapmt(): Unit {
    @Expect(dlapmtTest(lFloat64Implementation))
}

@Test
@Parallel
public func testDsytd2(): Unit {
    @Expect(dsytd2Test(lFloat64Implementation))
}

@Test
@Parallel
public func testDlangb(): Unit {
    @Expect(dlangbTest(lFloat64Implementation))
}

@Test
@Parallel
public func testDtrevc3(): Unit { 
    @Expect(dtrevc3Test(lFloat64Implementation))
}

@Test
@Parallel
public func testDpstrf(): Unit { // TODO ok?
    @Expect(dpstrfTest(lFloat64Implementation))
}

@Test
@Parallel
public func testDgelqf(): Unit { // TODO ok?
    @Expect(dgelqfTest(lFloat64Implementation))
}

@Test
@Parallel
public func testDgehrd(): Unit { // TODO ok?
    @Expect(dgehrdTest(lFloat64Implementation))
}

@Test
@Parallel
public func testDormr2(): Unit { // TODO ok?
    @Expect(dormr2Test(lFloat64Implementation))
}

@Test
@Parallel
public func testDgerqf(): Unit { // TODO ok ?
    @Expect(dgerqfTest(lFloat64Implementation))
}

@Test
@Parallel
public func testDorml2(): Unit { // TODO ok?
    @Expect(dorml2Test(lFloat64Implementation))
}

@Test
@Parallel
public func testDormlq(): Unit { // TODO ok?
    @Expect(dormlqTest(lFloat64Implementation))
}

@Test
@Parallel
public func testDormbr(): Unit { // TODO ok ?
    @Expect(dormbrTest(lFloat64Implementation))
}

@Test
@Skip // 精度 go 的精度也不够
public func testDpttrs(): Unit { // TODO
    @Expect(dpttrsTest(lFloat64Implementation))
}

@Test
@Parallel
public func testDlaqr1(): Unit {
    @Expect(dlaqr1Test(lFloat64Implementation))
}

@Test
@Parallel
public func testDtrtri(): Unit { 
    @Expect(dtrtriTest(lFloat64Implementation))
}

@Test
@Skip // 内存
public func testDgetrf(): Unit { // TODO
    @Expect(dgetrfTest(lFloat64Implementation))
}

@Test
@Parallel
public func testDgecon(): Unit { // TODO ok ?
    @Expect(dgeconTest(lFloat64Implementation))
}

@Test
@Parallel
public func testDgetri(): Unit { // TODO ok ?
    @Expect(dgetriTest(lFloat64Implementation))
}

@Test
@Parallel
public func testDgetrs(): Unit { // TODO ok ?
    @Expect(dgetrsTest(lFloat64Implementation))
}

@Test
@Parallel
public func testDlaqr04(): Unit { // TODO ok ?
    @Expect(dlaqr04Test(lFloat64Implementation))
}

@Test
@Parallel
public func testDlaqr23(): Unit { // TODO ok ?
    @Expect(dlaqr23Test(lFloat64Implementation))
}

@Test
@Parallel
public func testDormqr(): Unit {
    @Expect(dormqrTest(lFloat64Implementation))
}

@Test
@Parallel
public func testDorm2r(): Unit {
    @Expect(dorm2rTest(lFloat64Implementation))
}

@Test
@Parallel
public func testDormhr(): Unit { // TODO ok?
    @Expect(dormhrTest(lFloat64Implementation))
}

@Test
@Parallel
public func testDgeqrf(): Unit {
    @Expect(dgeqrfTest(lFloat64Implementation))
}

@Test
@Parallel
public func testDgebrd(): Unit { 
    @Expect(dgebrdTest(lFloat64Implementation))
}

@Test
@Parallel
public func testDlabrd(): Unit {
    @Expect(dlabrdTest(lFloat64Implementation))
}

@Test
@Skip // 依赖 dlasq2
public func testDlasq1(): Unit { // TODO ok ?
    @Expect(dlasq1Test(lFloat64Implementation))
}

@Test
@Skip // 精度
public func testDlasq2(): Unit { // TODO
    @Expect(dlasq2Test(lFloat64Implementation))
}

@Test
@Skip
public func testDorgql(): Unit { // TODO ok？
    @Expect(dorgqlTest(lFloat64Implementation))
}

@Test
@Skip
public func testDorgtr(): Unit { // TODO ok ?
    @Expect(dorgtrTest(lFloat64Implementation))
}

@Test
@Parallel
public func testDsteqr(): Unit { // TODO ok？
    @Expect(dsteqrTest(lFloat64Implementation))
}

@Test
@Parallel
public func testDsterf(): Unit { // TODO ok ?
    @Expect(dsterfTest(lFloat64Implementation))
}

@Test
@Parallel
public func testDsyev(): Unit { // TODO ok ?
    @Expect(dsyevTest(lFloat64Implementation))
}

@Test
@Parallel
public func testDsytrd(): Unit { // TODO ok ?
    @Expect(dsytrdTest(lFloat64Implementation))
}

@Test
@Parallel
public func testDorg2r(): Unit {
    @Expect(dorg2rTest(lFloat64Implementation))
}

@Test
@Parallel
public func testDorgqr(): Unit {
    @Expect(dorgqrTest(lFloat64Implementation))
}
