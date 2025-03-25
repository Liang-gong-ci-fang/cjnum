package cjnum.internal.asm.f64

import std.unittest.*
import std.unittest.testmacro.*
import cjnum.floats.scalar.*

class l1NormTestStruct {
    l1NormTestStruct(
        let want: Float64,
        let x: Array<Float64>
    ) {}
}

@Test
public func testL1Norm() {
    var src_gd: Float64 = 1.0
    let l1NormTests = [
        l1NormTestStruct(0.0, Array<Float64>([])),
        l1NormTestStruct(2.0, Array<Float64>([2.0])),
        l1NormTestStruct(6.0, Array<Float64>([1.0, 2.0, 3.0])),
        l1NormTestStruct(6.0, Array<Float64>([-1.0, -2.0, -3.0])),
        l1NormTestStruct(nan, Array<Float64>([nan])),
        l1NormTestStruct(40.0, Array<Float64>([8.0, -8.0, 8.0, -8.0, 8.0])),
        l1NormTestStruct(5.0, Array<Float64>([0.0, 1.0, 0.0, -1.0, 0.0, 1.0, 0.0, -1.0, 0.0, 1.0]))
    ]

    for (j in 0..l1NormTests.size) {
        let v = l1NormTests[j]
        let g_ln: Int64 = 4 + j % 2
        let guardedX = guardVector(v.x, src_gd, g_ln)
        let src = guardedX[g_ln..guardedX.size - g_ln]
        let ret = l1Norm(src)
        @Expect(same(ret, v.want))
        @Expect(isValidGuard(guardedX, src_gd, g_ln))
    }
}

class l1NormIncTestStruct {
    l1NormIncTestStruct(
        let inc: Int64,
        let want: Float64,
        let x: Array<Float64>
    ) {}
}

@Test
public func testL1NormInc() {
    var src_gd: Float64 = 1.0
    let l1NormIncTests = [
        l1NormIncTestStruct(2, 0.0, Array<Float64>([])),
        l1NormIncTestStruct(3, 2.0, Array<Float64>([2.0])),
        l1NormIncTestStruct(10, 6.0, Array<Float64>([1.0, 2.0, 3.0])),
        l1NormIncTestStruct(5, 6.0, Array<Float64>([-1.0, -2.0, -3.0])),
        l1NormIncTestStruct(3, nan, Array<Float64>([nan])),
        l1NormIncTestStruct(15, 40.0, Array<Float64>([8.0, -8.0, 8.0, -8.0, 8.0])),
        l1NormIncTestStruct(1, 5.0, Array<Float64>([0.0, 1.0, 0.0, -1.0, 0.0, 1.0, 0.0, -1.0, 0.0, 1.0]))
    ]

    for (j in 0..l1NormIncTests.size) {
        let v = l1NormIncTests[j]
        let g_ln: Int64 = 4 + j % 2
        let ln = v.x.size
        let guardedX = guardIncVector(v.x, src_gd, v.inc, g_ln)
        let src = guardedX[g_ln..guardedX.size - g_ln]
        let ret = l1NormInc(src, ln, v.inc)
        @Expect(same(ret, v.want))
        @Expect(checkValidIncGuard(guardedX, src_gd, v.inc, g_ln))
    }
}

class addTestStruct {
    addTestStruct(
        let dst: Array<Float64>,
        let src: Array<Float64>,
        let expect: Array<Float64>
    ) {}
}

@Test
public func testAdd() {
    let src_gd: Float64 = 1.0
    let dst_gd: Float64 = 0.0
    let addTests = [
        addTestStruct(Array<Float64>([1.0]), Array<Float64>([0.0]), Array<Float64>([1.0])),
        addTestStruct(Array<Float64>([1.0, 2.0, 3.0]), Array<Float64>([1.0]), Array<Float64>([2.0, 2.0, 3.0])),
        addTestStruct(Array<Float64>([]), Array<Float64>([]), Array<Float64>([])),
        addTestStruct(Array<Float64>([1.0]), Array<Float64>([nan]), Array<Float64>([nan])),
        addTestStruct(Array<Float64>([8.0, 8.0, 8.0, 8.0, 8.0]), Array<Float64>([2.0, 4.0, nan, 8.0, 9.0]),
            Array<Float64>([10.0, 12.0, nan, 16.0, 17.0])),
        addTestStruct(Array<Float64>([0.0, 1.0, 2.0, 3.0, 4.0]), Array<Float64>([-inf, 4.0, nan, 8.0, 9.0]),
            Array<Float64>([-inf, 5.0, nan, 11.0, 13.0])),
        addTestStruct(Array<Float64>(48, item: 0.0), Array<Float64>(48, item: 0.0), Array<Float64>(48, item: 0.0))
    ]

    for (j in 0..addTests.size) {
        let v = addTests[j]
        let sg_ln: Int64 = 4 + j % 2
        let dg_ln: Int64 = 4 + j % 3
        let guardedSrc = guardVector(v.src, src_gd, sg_ln)
        let guardedDst = guardVector(v.dst, dst_gd, dg_ln)
        let src = guardedSrc[sg_ln..guardedSrc.size - sg_ln]
        let dst = guardedDst[dg_ln..guardedDst.size - dg_ln]
        add(dst, src)
        for (i in 0..v.expect.size) {
            @Expect(same(dst[i], v.expect[i]))
        }
        @Expect(isValidGuard(guardedSrc, src_gd, sg_ln))
        @Expect(isValidGuard(guardedDst, dst_gd, dg_ln))
    }
}

class addConstTestStruct {
    addConstTestStruct(
        let alpha: Float64,
        let src: Array<Float64>,
        let expect: Array<Float64>
    ) {}
}

@Test
public func testAddConst() {
    let src_gd: Float64 = 0.0
    let addConstTests = [
        addConstTestStruct(1.0, Array<Float64>([0.0]), Array<Float64>([1.0])),
        addConstTestStruct(5.0, Array<Float64>([]), Array<Float64>([])),
        addConstTestStruct(1.0, Array<Float64>([nan]), Array<Float64>([nan])),
        addConstTestStruct(8.0, Array<Float64>()),
        addConstTestStruct(inf, Array<Float64>([-inf, 4.0, nan, 8.0, 9.0]), Array<Float64>([nan, inf, nan, inf, inf]))
    ]

    for (j in 0..addConstTests.size) {
        let v = addConstTests[j]
        let g_ln: Int64 = 4 + j % 2
        let guardedSrc = guardVector(v.src, src_gd, g_ln)
        let src = guardedSrc[g_ln..guardedSrc.size - g_ln]
        addConst(v.alpha, src)
        for (i in 0..v.expect.size) {
            @Expect(same(src[i], v.expect[i]))
        }
        @Expect(isValidGuard(guardedSrc, src_gd, g_ln))
    }
}

class cumSumTestStruct {
    cumSumTestStruct(
        let dst: Array<Float64>,
        let src: Array<Float64>,
        let expect: Array<Float64>
    ) {}
}

@Test
public func testCumSum() {
    let src_gd: Float64 = -1.0
    let dst_gd: Float64 = 0.0
    let cumSumTests = [
        cumSumTestStruct(Array<Float64>([]), Array<Float64>([]), Array<Float64>([])),
        cumSumTestStruct(Array<Float64>([0.0]), Array<Float64>([1.0]), Array<Float64>([1.0])),
        cumSumTestStruct(Array<Float64>([nan]), Array<Float64>([nan]), Array<Float64>([nan])),
        cumSumTestStruct(Array<Float64>([0.0, 0.0, 0.0]), Array<Float64>([1.0, 2.0, 3.0]), Array<Float64>([1.0, 3.0, 6.0]
        )),
        cumSumTestStruct(Array<Float64>([0.0, 0.0, 0.0, 0.0]), Array<Float64>([1.0, 2.0, 3.0]),
            Array<Float64>([1.0, 3.0, 6.0])),
        cumSumTestStruct(Array<Float64>([0.0, 0.0, 0.0, 0.0]), Array<Float64>([1.0, 2.0, 3.0, 4.0]),
            Array<Float64>([1.0, 3.0, 6.0, 10.0])),
        cumSumTestStruct(Array<Float64>([1.0, nan, nan, 1.0, 1.0]), Array<Float64>([1.0, 1.0, nan, 1.0, 1.0]),
            Array<Float64>([1.0, 2.0, nan, nan, nan])),
        cumSumTestStruct(Array<Float64>([nan, 4.0, inf, -inf, 9.0]), Array<Float64>([inf, 4.0, nan, -inf, 9.0]),
            Array<Float64>([inf, inf, nan, nan, nan])),
        cumSumTestStruct(Array<Float64>(16, item: 0.0),
            Array<Float64>([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]),
            Array<Float64>([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0]))
    ]

    for (j in 0..cumSumTests.size) {
        let v = cumSumTests[j]
        let g_ln: Int64 = 4 + j % 2
        let guardedSrc = guardVector(v.src, src_gd, g_ln)
        let guardedDst = guardVector(v.dst, dst_gd, g_ln)
        let src = guardedSrc[g_ln..guardedSrc.size - g_ln]
        let dst = guardedDst[g_ln..guardedDst.size - g_ln]
        cumSum(dst, src)
        for (i in 0..v.expect.size) {
            @Expect(same(dst[i], v.expect[i]))
        }
        @Expect(isValidGuard(guardedSrc, src_gd, g_ln))
        @Expect(isValidGuard(guardedDst, dst_gd, g_ln))
    }
}

class cumProdTestStruct {
    cumProdTestStruct(
        let dst: Array<Float64>,
        let src: Array<Float64>,
        let expect: Array<Float64>
    ) {}
}

@Test
public func testCumProd() {
    let src_gd: Float64 = -1.0
    let dst_gd: Float64 = 1.0
    let cumProdTests = [
        cumProdTestStruct(Array<Float64>([]), Array<Float64>([]), Array<Float64>([])),
        cumProdTestStruct(Array<Float64>([1.0]), Array<Float64>([1.0]), Array<Float64>([1.0])),
        cumProdTestStruct(Array<Float64>([nan]), Array<Float64>([nan]), Array<Float64>([nan])),
        cumProdTestStruct(Array<Float64>([0.0, 0.0, 0.0, 0.0]), Array<Float64>([1.0, 2.0, 3.0, 4.0]),
            Array<Float64>([1.0, 2.0, 6.0, 24.0])),
        cumProdTestStruct(Array<Float64>([0.0, 0.0, 0.0]), Array<Float64>([1.0, 2.0, 3.0]),
            Array<Float64>([1.0, 2.0, 6.0])),
        cumProdTestStruct(Array<Float64>([0.0, 0.0, 0.0, 0.0]), Array<Float64>([1.0, 2.0, 3.0]),
            Array<Float64>([1.0, 2.0, 6.0])),
        cumProdTestStruct(Array<Float64>([nan, 1.0, nan, 1.0, 0.0]), Array<Float64>([1.0, 1.0, nan, 1.0, 1.0]),
            Array<Float64>([1.0, 1.0, nan, nan, nan])),
        cumProdTestStruct(Array<Float64>([nan, 4.0, nan, -inf, 9.0]), Array<Float64>([inf, 4.0, nan, -inf, 9.0]),
            Array<Float64>([inf, inf, nan, nan, nan])),
        cumProdTestStruct(Array<Float64>(18, item: 0.0),
            Array<Float64>([2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0]),
            Array<Float64>(
            [2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0, 256.0, 512.0, 1024.0, 2048.0, 4096.0, 8192.0, 16384.0, 32768.0,
            65536.0]))
    ]

    for (j in 0..cumProdTests.size) {
        let v = cumProdTests[j]
        let sg_ln: Int64 = 4 + j % 2
        let dg_ln: Int64 = 4 + j % 3
        let guardedSrc = guardVector(v.src, src_gd, sg_ln)
        let guardedDst = guardVector(v.dst, dst_gd, dg_ln)
        let src = guardedSrc[sg_ln..guardedSrc.size - sg_ln]
        let dst = guardedDst[dg_ln..guardedDst.size - dg_ln]
        cumProd(dst, src)
        for (i in 0..v.expect.size) {
            @Expect(same(dst[i], v.expect[i]))
        }
        @Expect(isValidGuard(guardedSrc, src_gd, sg_ln))
        @Expect(isValidGuard(guardedDst, dst_gd, dg_ln))
    }
}

class divTestStruct {
    divTestStruct(
        let dst: Array<Float64>,
        let src: Array<Float64>,
        let expect: Array<Float64>
    ) {}
}

@Test
public func testDiv() {
    let src_gd: Float64 = -1.0
    let dst_gd: Float64 = 0.5
    let divTests = [
        divTestStruct(Array<Float64>([1.0]), Array<Float64>([1.0]), Array<Float64>([1.0])),
        divTestStruct(Array<Float64>([nan]), Array<Float64>([nan]), Array<Float64>([nan])),
        divTestStruct(Array<Float64>([1.0, 2.0, 3.0, 4.0]), Array<Float64>([1.0, 2.0, 3.0, 4.0]),
            Array<Float64>([1.0, 1.0, 1.0, 1.0])),
        divTestStruct(Array<Float64>([1.0, 2.0, 3.0, 4.0, 2.0, 4.0, 6.0, 8.0]),
            Array<Float64>([1.0, 2.0, 3.0, 4.0, 1.0, 2.0, 3.0, 4.0]), Array<Float64>([1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0,
            2.0])),
        divTestStruct(Array<Float64>([2.0, 4.0, 6.0]), Array<Float64>([1.0, 2.0, 3.0]), Array<Float64>([2.0, 2.0, 2.0])),
        divTestStruct(Array<Float64>([0.0, 0.0, 0.0, 0.0]), Array<Float64>([1.0, 2.0, 3.0]),
            Array<Float64>([0.0, 0.0, 0.0])),
        divTestStruct(Array<Float64>([nan, 1.0, nan, 1.0, 0.0, nan, 1.0, nan, 1.0, 0.0]),
            Array<Float64>([1.0, 1.0, nan, 1.0, 1.0, 1.0, 1.0, nan, 1.0, 1.0]),
            Array<Float64>([nan, 1.0, nan, 1.0, 0.0, nan, 1.0, nan, 1.0, 0.0])),
        divTestStruct(Array<Float64>([inf, 4.0, nan, -inf, 9.0, inf, 4.0, nan, -inf, 9.0]),
            Array<Float64>([inf, 4.0, nan, -inf, 3.0, inf, 4.0, nan, -inf, 3.0]),
            Array<Float64>([nan, 1.0, nan, nan, 3.0, nan, 1.0, nan, nan, 3.0]))
    ]

    for (j in 0..divTests.size) {
        let v = divTests[j]
        let sg_ln: Int64 = 4 + j % 2
        let dg_ln: Int64 = 4 + j % 3
        let guardedSrc = guardVector(v.src, src_gd, sg_ln)
        let guardedDst = guardVector(v.dst, dst_gd, dg_ln)
        let src = guardedSrc[sg_ln..guardedSrc.size - sg_ln]
        let dst = guardedDst[dg_ln..guardedDst.size - dg_ln]
        div(dst, src)
        for (i in 0..v.expect.size) {
            @Expect(same(dst[i], v.expect[i]))
        }
        @Expect(isValidGuard(guardedSrc, src_gd, sg_ln))
        @Expect(isValidGuard(guardedDst, dst_gd, dg_ln))
    }
}

class divToTestStruct {
    divToTestStruct(
        let dst: Array<Float64>,
        let x: Array<Float64>,
        let y: Array<Float64>,
        let expect: Array<Float64>
    ) {}
}

@Test
public func testDivTo() {
    let dst_gd: Float64 = -1.0
    let x_gd: Float64 = 0.5
    let y_gd: Float64 = 0.25
    let divToTests = [
        divToTestStruct(Array<Float64>([1.0]), Array<Float64>([1.0]), Array<Float64>([1.0]), Array<Float64>([1.0])),
        divToTestStruct(Array<Float64>([1.0]), Array<Float64>([nan]), Array<Float64>([nan]), Array<Float64>([nan])),
        divToTestStruct(Array<Float64>([-2.0, -2.0, -2.0]), Array<Float64>([1.0, 2.0, 3.0]),
            Array<Float64>([1.0, 2.0, 3.0]), Array<Float64>([1.0, 1.0, 1.0])),
        divToTestStruct(Array<Float64>([0.0, 0.0, 0.0]), Array<Float64>([2.0, 4.0, 6.0]),
            Array<Float64>([1.0, 2.0, 3.0, 4.0]), Array<Float64>([2.0, 2.0, 2.0])),
        divToTestStruct(Array<Float64>([-1.0, -1.0, -1.0]), Array<Float64>([0.0, 0.0, 0.0]),
            Array<Float64>([1.0, 2.0, 3.0]), Array<Float64>([0.0, 0.0, 0.0])),
        divToTestStruct(Array<Float64>([inf, inf, inf, inf, inf, inf, inf, inf, inf, inf]),
            Array<Float64>([nan, 1.0, nan, 1.0, 0.0, nan, 1.0, nan, 1.0, 0.0]),
            Array<Float64>([1.0, 1.0, nan, 1.0, 1.0, 1.0, 1.0, nan, 1.0, 1.0]),
            Array<Float64>([nan, 1.0, nan, 1.0, 0.0, nan, 1.0, nan, 1.0, 0.0])),
        divToTestStruct(Array<Float64>([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
            Array<Float64>([inf, 4.0, nan, -inf, 9.0, inf, 4.0, nan, -inf, 9.0]),
            Array<Float64>([inf, 4.0, nan, -inf, 3.0, inf, 4.0, nan, -inf, 3.0]),
            Array<Float64>([nan, 1.0, nan, nan, 3.0, nan, 1.0, nan, nan, 3.0]))
    ]

    for (j in 0..divToTests.size) {
        let v = divToTests[j]
        let xg_ln: Int64 = 4 + j % 2
        let yg_ln: Int64 = 4 + j % 3
        let guardedX = guardVector(v.x, x_gd, xg_ln)
        let guardedY = guardVector(v.y, y_gd, yg_ln)
        let x = guardedX[xg_ln..guardedX.size - xg_ln]
        let y = guardedY[yg_ln..guardedY.size - yg_ln]
        let guardedDst = guardVector(v.dst, dst_gd, xg_ln)
        let dst = guardedDst[xg_ln..guardedDst.size - xg_ln]
        let ret = divTo(dst, x, y)
        for (i in 0..v.expect.size) {
            @Expect(same(ret[i], v.expect[i]))
            @Expect(same(ret[i], dst[i]))
        }
        @Expect(isValidGuard(guardedX, x_gd, xg_ln))
        @Expect(isValidGuard(guardedY, y_gd, yg_ln))
        @Expect(isValidGuard(guardedDst, dst_gd, xg_ln))
    }
}

class L1DistTestStruct {
    L1DistTestStruct(
        let s: Array<Float64>,
        let t: Array<Float64>,
        let expect: Float64
    ) {}
}

@Test
public func testL1Dist() {
    let s_gd: Float64 = -1.0
    let t_gd: Float64 = 0.5
    let l1DistTests = [
        L1DistTestStruct(Array<Float64>([1.0]), Array<Float64>([1.0]), 0.0),
        L1DistTestStruct(Array<Float64>([nan]), Array<Float64>([nan]), nan),
        L1DistTestStruct(Array<Float64>([1.0, 2.0, 3.0, 4.0]), Array<Float64>([1.0, 2.0, 3.0, 4.0]), 0.0),
        L1DistTestStruct(Array<Float64>([2.0, 4.0, 6.0]), Array<Float64>([1.0, 2.0, 3.0, 4.0]), 6.0),
        L1DistTestStruct(Array<Float64>([0.0, 0.0, 0.0]), Array<Float64>([1.0, 2.0, 3.0]), 6.0),
        L1DistTestStruct(Array<Float64>([0.0, -4.0, -10.0]), Array<Float64>([1.0, 2.0, 3.0]), 20.0),
        L1DistTestStruct(Array<Float64>([0.0, 1.0, 0.0, 1.0, 0.0]), Array<Float64>([1.0, 1.0, inf, 1.0, 1.0]), inf),
        L1DistTestStruct(Array<Float64>([inf, 4.0, nan, -inf, 9.0]), Array<Float64>([inf, 4.0, nan, -inf, 3.0]), nan)
    ]

    for (j in 0..l1DistTests.size) {
        let v = l1DistTests[j]
        let sg_ln: Int64 = 4 + j % 2
        let tg_ln: Int64 = 4 + j % 3
        let guardedS = guardVector(v.s, s_gd, sg_ln)
        let guardedT = guardVector(v.t, t_gd, tg_ln)
        let s = guardedS[sg_ln..guardedS.size - sg_ln]
        let t = guardedT[tg_ln..guardedT.size - tg_ln]
        let ret = l1Dist(s, t)
        @Expect(same(ret, v.expect))
        @Expect(isValidGuard(guardedS, s_gd, sg_ln))
        @Expect(isValidGuard(guardedT, t_gd, tg_ln))
    }
}

class linfDistTestStruct {
    linfDistTestStruct(
        let s: Array<Float64>,
        let t: Array<Float64>,
        let expect: Float64
    ) {}
}

@Test
public func testLinfDist() {
    let t_gd: Float64 = 0.0
    let s_gd: Float64 = inf
    let linfDistTests = [
        linfDistTestStruct(Array<Float64>([]), Array<Float64>([]), 0.0),
        linfDistTestStruct(Array<Float64>([1.0]), Array<Float64>([1.0]), 0.0),
        linfDistTestStruct(Array<Float64>([nan]), Array<Float64>([nan]), nan),
        linfDistTestStruct(Array<Float64>([1.0, 2.0, 3.0, 4.0]), Array<Float64>([1.0, 2.0, 3.0, 4.0]), 0.0),
        linfDistTestStruct(Array<Float64>([2.0, 4.0, 6.0]), Array<Float64>([1.0, 2.0, 3.0, 4.0]), 3.0),
        linfDistTestStruct(Array<Float64>([0.0, 0.0, 0.0]), Array<Float64>([1.0, 2.0, 3.0]), 3.0),
        linfDistTestStruct(Array<Float64>([0.0, 1.0, 0.0, 1.0, 0.0]), Array<Float64>([1.0, 1.0, inf, 1.0, 1.0]), inf),
        linfDistTestStruct(Array<Float64>([inf, 4.0, nan, -inf, 9.0]), Array<Float64>([inf, 4.0, nan, -inf, 3.0]), 6.0)
    ]

    for (j in 0..linfDistTests.size) {
        let v = linfDistTests[j]
        let sg_ln: Int64 = 4 + j % 2
        let tg_ln: Int64 = 4 + j % 3
        let guardedS = guardVector(v.s, s_gd, sg_ln)
        let guardedT = guardVector(v.t, t_gd, tg_ln)
        let s = guardedS[sg_ln..guardedS.size - sg_ln]
        let t = guardedT[tg_ln..guardedT.size - tg_ln]
        let ret = linfDist(s, t)
        @Expect(same(ret, v.expect))
        @Expect(isValidGuard(guardedS, s_gd, sg_ln))
        @Expect(isValidGuard(guardedT, t_gd, tg_ln))
    }
}

class SumTestStruct {
    SumTestStruct(
        let src: Array<Float64>,
        let expect: Float64
    ) {}
}

@Test
public func testSum() {
    let src_gd: Float64 = -1.0
    let sumTests = [
        SumTestStruct(Array<Float64>([]), 0.0),
        SumTestStruct(Array<Float64>([1.0]), 1.0),
        SumTestStruct(Array<Float64>([nan]), nan),
        SumTestStruct(Array<Float64>([1.0, 2.0, 3.0]), 6.0),
        SumTestStruct(Array<Float64>([1.0, -4.0, 3.0]), 0.0),
        SumTestStruct(Array<Float64>([1.0, 2.0, 3.0, 4.0]), 10.0),
        SumTestStruct(Array<Float64>([1.0, 1.0, nan, 1.0, 1.0]), nan),
        SumTestStruct(Array<Float64>([inf, 4.0, nan, -inf, 9.0]), nan),
        SumTestStruct(Array<Float64>([1.0, 1.0, 1.0, 1.0, 9.0, 1.0, 1.0, 1.0, 2.0, 1.0, 1.0, 1.0, 1.0, 1.0, 5.0, 1.0]),
            29.0),
        SumTestStruct(
            Array<Float64>(
            [1.0, 1.0, 1.0, 1.0, 9.0, 1.0, 1.0, 1.0, 2.0, 1.0, 1.0, 1.0, 1.0, 1.0, 5.0, 11.0, 1.0, 1.0, 1.0, 9.0, 1.0,
            1.0, 1.0, 2.0, 1.0, 1.0, 1.0, 1.0, 1.0, 5.0, 1.0]), 67.0)
    ]

    for (j in 0..sumTests.size) {
        let v = sumTests[j]
        var gdLn: Int64 = 4 + j % 2
        var gsrc = guardVector(v.src, src_gd, gdLn)
        var src = gsrc[gdLn..gsrc.size - gdLn]
        var ret = sum(src)
        @Expect(same(ret, v.expect))
        @Expect(isValidGuard(gsrc, src_gd, gdLn))

        gdLn++
        gsrc = guardVector(v.src, src_gd, gdLn)
        src = gsrc[gdLn..gsrc.size - gdLn]
        ret = sum(src)
        @Expect(same(ret, v.expect))
        @Expect(isValidGuard(gsrc, src_gd, gdLn))
    }
}
