package cjnum.internal.asm.c128

import cjnum.cmplxs.cscalar.*
import cjnum.complex.*
import std.math.*
import std.unittest.*
import std.unittest.testmacro.*

class addTestStruct {
    addTestStruct(
        let dst!: Array<Complex128>,
        let src!: Array<Complex128>,
        let expect!: Array<Complex128>
    ) {}
}

@Test
public func testAdd() {
    var src_gd: Complex128 = Complex128(1.0)
    var dst_gd: Complex128 = Complex128(0.0)
    let tests = [
        addTestStruct(dst: [Complex128(1.0, 1.0)], src: [Complex128(0.0)], expect: [Complex128(1.0, 1.0)]),
        addTestStruct(
            dst: [Complex128(1.0, 0.0), Complex128(2.0, 0.0), Complex128(3.0, 0.0)],
            src: [Complex128(1.0, 1.0)],
            expect: [Complex128(2.0, 1.0), Complex128(2.0, 0.0), Complex128(3.0, 0.0)]
        ),
        addTestStruct(
            dst: [],
            src: [],
            expect: []
        ),
        addTestStruct(
            dst: [Complex128(1.0, 0.0)],
            src: [Complex128.NaN],
            expect: [Complex128.NaN]
        ),
        addTestStruct(
            dst: [Complex128(8.0, 0.0), Complex128(8.0, 0.0), Complex128(8.0, 0.0), Complex128(8.0, 0.0),
                Complex128(8.0, 0.0)],
            src: [Complex128(2.0, 1.0), Complex128(4.0, -1.0), Complex128.NaN, Complex128(8.0, 1.0),
                Complex128(9.0, -1.0)],
            expect: [Complex128(10.0, 1.0), Complex128(12.0, -1.0), Complex128.NaN, Complex128(16.0, 1.0),
                Complex128(17.0, -1.0)]
        ),
        addTestStruct(
            dst: [Complex128(0.0, 0.0), Complex128(1.0, 1.0), Complex128(2.0, 0.0), Complex128(3.0, -1.0),
                Complex128(4.0, 0.0)],
            src: [Complex128(inf), Complex128(4.0, 0.0), Complex128.NaN, Complex128(8.0, 1.0), Complex128(9.0, -1.0)],
            expect: [Complex128(inf), Complex128(5.0, 1.0), Complex128.NaN, Complex128(11.0, 0.0), Complex128(13.0, -1.0
            )]
        ),
        addTestStruct(
            dst: Array<Complex128>(48, item: Complex128(0.0, 0.0)),
            src: Array<Complex128>(48, item: Complex128(0.0, 0.0)),
            expect: Array<Complex128>(48, item: Complex128(0.0, 0.0))
        )
    ]
    for (j in 0..tests.size) {
        let sg_ln = 4 + j % 2
        let dg_ln = 4 + j % 3
        let src_ = guardVector(tests[j].src, src_gd, sg_ln)
        let dst_ = guardVector(tests[j].dst, dst_gd, dg_ln)
        let src = src_[sg_ln..src_.size - sg_ln]
        let dst = dst_[dg_ln..dst_.size - dg_ln]
        add(dst, src)
        for (i in 0..tests[j].expect.size) {
            @Expect(same(dst[i], tests[j].expect[i]))
        }
        @Expect(isValidGuard(src_, src_gd, sg_ln))
        @Expect(isValidGuard(dst_, dst_gd, dg_ln))
    }
}

class addConstTestStruct {
    addConstTestStruct(
        let alpha!: Complex128,
        let src!: Array<Complex128>,
        let expect!: Array<Complex128>
    ) {}
}

@Test
public func testAddConst() {
    var src_gd: Complex128 = Complex128(0.0)
    let tests = [
        addConstTestStruct(alpha: Complex128(1.0, 1.0), src: [Complex128(0.0)], expect: [Complex128(1.0, 1.0)]),
        addConstTestStruct(alpha: Complex128(5.0), src: [], expect: []),
        addConstTestStruct(alpha: Complex128(1.0), src: [Complex128.NaN], expect: [Complex128.NaN]),
        addConstTestStruct(
            alpha: Complex128(8.0, 1.0),
            src: [Complex128(2.0), Complex128(4.0), Complex128.NaN, Complex128(8.0), Complex128(9.0)],
            expect: [Complex128(10.0, 1.0), Complex128(12.0, 1.0), Complex128.NaN, Complex128(16.0, 1.0),
                Complex128(17.0, 1.0)]
        ),
        addConstTestStruct(
            alpha: Complex128(inf),
            src: [Complex128(inf), Complex128(4.0), Complex128.NaN, Complex128(8.0), Complex128(9.0)],
            expect: [Complex128(inf), Complex128(inf), Complex128.NaN, Complex128(inf), Complex128(inf)]
        )
    ]
    for (j in 0..tests.size) {
        let g_ln = 4 + j % 2
        let src_ = guardVector(tests[j].src, src_gd, g_ln)
        let src = src_[g_ln..src_.size - g_ln]
        addConst(tests[j].alpha, src)
        for (i in 0..tests[j].expect.size) {
            @Expect(same(src[i], tests[j].expect[i]))
        }
        @Expect(isValidGuard(src_, src_gd, g_ln))
    }
}

class cumSumTestStruct {
    cumSumTestStruct(
        let dst!: Array<Complex128>,
        let src!: Array<Complex128>,
        let expect!: Array<Complex128>
    ) {}
}

@Test
public func testCumSum() {
    var src_gd = Complex128(-1.0)
    var dst_gd = Complex128(0.0)
    let tests = [
        cumSumTestStruct(dst: [], src: [], expect: []),
        cumSumTestStruct(dst: [Complex128(0.0)], src: [Complex128(1.0, 1.0)], expect: [Complex128(1.0, 1.0)]),
        cumSumTestStruct(dst: [Complex128.NaN], src: [Complex128.NaN], expect: [Complex128.NaN]),
        cumSumTestStruct(
            dst: [Complex128(0.0), Complex128(0.0), Complex128(0.0)],
            src: [Complex128(1.0), Complex128(2.0, 1.0), Complex128(3.0, 2.0)],
            expect: [Complex128(1.0), Complex128(3.0, 1.0), Complex128(6.0, 3.0)]
        ),
        cumSumTestStruct(
            dst: [Complex128(0.0), Complex128(0.0), Complex128(0.0), Complex128(0.0)],
            src: [Complex128(1.0), Complex128(2.0, 1.0), Complex128(3.0, 2.0)],
            expect: [Complex128(1.0), Complex128(3.0, 1.0), Complex128(6.0, 3.0)]
        ),
        cumSumTestStruct(
            dst: [Complex128(0.0), Complex128(0.0), Complex128(0.0), Complex128(0.0)],
            src: [Complex128(1.0), Complex128(2.0, 1.0), Complex128(3.0, 2.0), Complex128(4.0, 3.0)],
            expect: [Complex128(1.0), Complex128(3.0, 1.0), Complex128(6.0, 3.0), Complex128(10.0, 6.0)]
        ),
        cumSumTestStruct(
            dst: [Complex128(1.0), Complex128.NaN, Complex128.NaN, Complex128(1.0), Complex128(1.0)],
            src: [Complex128(1.0), Complex128(1.0), Complex128.NaN, Complex128(1.0), Complex128(1.0)],
            expect: [Complex128(1.0), Complex128(2.0), Complex128.NaN, Complex128.NaN, Complex128.NaN]
        ),
        cumSumTestStruct(
            dst: [Complex128.NaN, Complex128(4.0), Complex128(inf), Complex128(inf), Complex128(9.0)],
            src: [Complex128(inf), Complex128(4.0), Complex128.NaN, Complex128(inf), Complex128(9.0)],
            expect: [Complex128(inf), Complex128(inf), Complex128.NaN, Complex128.NaN, Complex128.NaN]
        ),
        cumSumTestStruct(
            dst: [Complex128(0.0), Complex128(0.0), Complex128(0.0), Complex128(0.0), Complex128(0.0), Complex128(0.0),
                Complex128(0.0), Complex128(0.0), Complex128(0.0), Complex128(0.0), Complex128(0.0), Complex128(0.0),
                Complex128(0.0), Complex128(0.0), Complex128(0.0), Complex128(0.0)],
            src: [Complex128(1.0), Complex128(1.0), Complex128(1.0), Complex128(1.0), Complex128(1.0), Complex128(1.0),
                Complex128(1.0), Complex128(1.0), Complex128(1.0), Complex128(1.0), Complex128(1.0), Complex128(1.0),
                Complex128(1.0), Complex128(1.0), Complex128(1.0), Complex128(1.0)],
            expect: [Complex128(1.0), Complex128(2.0), Complex128(3.0), Complex128(4.0), Complex128(5.0),
                Complex128(6.0), Complex128(7.0), Complex128(8.0), Complex128(9.0), Complex128(10.0), Complex128(11.0),
                Complex128(12.0), Complex128(13.0), Complex128(14.0), Complex128(15.0), Complex128(16.0)]
        )
    ]
    for (j in 0..tests.size) {
        let g_ln = 4 + j % 2
        let src_ = guardVector(tests[j].src, src_gd, g_ln)
        let dst_ = guardVector(tests[j].dst, dst_gd, g_ln)
        let src = src_[g_ln..src_.size - g_ln]
        let dst = dst_[g_ln..dst_.size - g_ln]
        let ret = cumSum(dst, src)
        for (i in 0..tests[j].expect.size) {
            @Expect(same(ret[i], tests[j].expect[i]))
            @Expect(same(ret[i], dst[i]))
        }
        @Expect(isValidGuard(src_, src_gd, g_ln))
        @Expect(isValidGuard(dst_, dst_gd, g_ln))
    }
}

class cumProdTestStruct {
    cumProdTestStruct(
        let dst!: Array<Complex128>,
        let src!: Array<Complex128>,
        let expect!: Array<Complex128>
    ) {}
}

@Test
public func testCumProd() {
    var src_gd = Complex128(-1.0)
    var dst_gd = Complex128(1.0)
    let tests = [
        cumProdTestStruct(dst: [], src: [], expect: []),
        cumProdTestStruct(dst: [Complex128(1.0)], src: [Complex128(1.0, 1.0)], expect: [Complex128(1.0, 1.0)]),
        cumProdTestStruct(dst: [Complex128.NaN], src: [Complex128.NaN], expect: [Complex128.NaN]),
        cumProdTestStruct(
            dst: [Complex128(0.0), Complex128(0.0), Complex128(0.0), Complex128(0.0)],
            src: [Complex128(1.0), Complex128(2.0, 1.0), Complex128(3.0, 2.0), Complex128(4.0, 3.0)],
            expect: [Complex128(1.0), Complex128(2.0, 1.0), Complex128(4.0, 7.0), Complex128(-5.0, 40.0)]
        ),
        cumProdTestStruct(
            dst: [Complex128(0.0), Complex128(0.0), Complex128(0.0)],
            src: [Complex128(1.0), Complex128(2.0, 1.0), Complex128(3.0, 2.0)],
            expect: [Complex128(1.0), Complex128(2.0, 1.0), Complex128(4.0, 7.0)]
        ),
        cumProdTestStruct(
            dst: [Complex128(0.0), Complex128(0.0), Complex128(0.0), Complex128(0.0)],
            src: [Complex128(1.0), Complex128(2.0, 1.0), Complex128(3.0, 2.0)],
            expect: [Complex128(1.0), Complex128(2.0, 1.0), Complex128(4.0, 7.0)]
        ),
        cumProdTestStruct(
            dst: [Complex128.NaN, Complex128(1.0), Complex128.NaN, Complex128(1.0), Complex128(0.0)],
            src: [Complex128(1.0), Complex128(1.0), Complex128.NaN, Complex128(1.0), Complex128(1.0)],
            expect: [Complex128(1.0), Complex128(1.0), Complex128.NaN, Complex128.NaN, Complex128.NaN]
        ),
        cumProdTestStruct(
            dst: Array<Complex128>(18, item: Complex128(0.0, 0.0)),
            src: [Complex128(0.0, 2.0), Complex128(0.0, 2.0), Complex128(0.0, 2.0), Complex128(0.0, 2.0),
                Complex128(0.0, 2.0), Complex128(0.0, 2.0), Complex128(0.0, 2.0), Complex128(0.0, 2.0),
                Complex128(0.0, 2.0), Complex128(0.0, 2.0), Complex128(0.0, 2.0), Complex128(0.0, 2.0),
                Complex128(0.0, 2.0), Complex128(0.0, 2.0), Complex128(0.0, 2.0), Complex128(0.0, 2.0),
                Complex128(0.0, 2.0), Complex128(0.0, 2.0)],
            expect: [Complex128(0.0, 2.0), Complex128(-4.0, 0.0), Complex128(0.0, -8.0), Complex128(16.0, 0.0),
                Complex128(0.0, 32.0), Complex128(-64.0, 0.0), Complex128(0.0, -128.0), Complex128(256.0, 0.0),
                Complex128(0.0, 512.0), Complex128(-1024.0, 0.0), Complex128(0.0, -2048.0), Complex128(4096.0, 0.0),
                Complex128(0.0, 8192.0), Complex128(-16384.0, 0.0), Complex128(0.0, -32768.0), Complex128(65536.0, 0.0)]
        )
    ]
    for (j in 0..tests.size) {
        let sg_ln = 4 + j % 2
        let dg_ln = 4 + j % 3
        let src_ = guardVector(tests[j].src, src_gd, sg_ln)
        let dst_ = guardVector(tests[j].dst, dst_gd, dg_ln)
        let src = src_[sg_ln..src_.size - sg_ln]
        let dst = dst_[dg_ln..dst_.size - dg_ln]
        let ret = cumProd(dst, src)
        for (i in 0..tests[j].expect.size) {
            @Expect(same(ret[i], tests[j].expect[i]))
            @Expect(same(ret[i], dst[i]))
        }
        @Expect(isValidGuard(src_, src_gd, sg_ln))
        @Expect(isValidGuard(dst_, dst_gd, dg_ln))
    }
}

class divTestStruct {
    divTestStruct(
        let dst!: Array<Complex128>,
        let src!: Array<Complex128>,
        let expect!: Array<Complex128>
    ) {}
}

@Test
public func testDiv() {
    const tol = 1e-15
    var src_gd = Complex128(-1.0)
    var dst_gd = Complex128(0.5)
    let tests = [
        divTestStruct(
            dst: [Complex128(1.0, 1.0)],
            src: [Complex128(1.0, 1.0)],
            expect: [Complex128(1.0)]
        ),
        divTestStruct(
            dst: [Complex128.NaN],
            src: [Complex128.NaN],
            expect: [Complex128.NaN]
        ),
        divTestStruct(
            dst: [Complex128(1.0, 1.0), Complex128(2.0, 1.0), Complex128(3.0, 1.0), Complex128(4.0, 1.0)],
            src: [Complex128(1.0, 1.0), Complex128(2.0, 1.0), Complex128(3.0, 1.0), Complex128(4.0, 1.0)],
            expect: [Complex128(1.0), Complex128(1.0), Complex128(1.0), Complex128(1.0)]
        ),
        divTestStruct(
            dst: [Complex128(1.0, 1.0), Complex128(2.0, 1.0), Complex128(3.0, 1.0), Complex128(4.0, 1.0),
                Complex128(2.0, 2.0), Complex128(4.0, 2.0), Complex128(6.0, 2.0), Complex128(8.0, 2.0)],
            src: [Complex128(1.0, 1.0), Complex128(2.0, 1.0), Complex128(3.0, 1.0), Complex128(4.0, 1.0),
                Complex128(1.0, 1.0), Complex128(2.0, 1.0), Complex128(3.0, 1.0), Complex128(4.0, 1.0)],
            expect: [Complex128(1.0), Complex128(1.0), Complex128(1.0), Complex128(1.0), Complex128(2.0),
                Complex128(2.0), Complex128(2.0), Complex128(2.0)]
        ),
        divTestStruct(
            dst: [Complex128(2.0, 2.0), Complex128(4.0, 8.0), Complex128(6.0, -12.0)],
            src: [Complex128(1.0, 1.0), Complex128(2.0, 4.0), Complex128(3.0, -6.0)],
            expect: [Complex128(2.0), Complex128(2.0), Complex128(2.0)]
        ),
        divTestStruct(
            dst: [Complex128(0.0), Complex128(0.0), Complex128(0.0), Complex128(0.0)],
            src: [Complex128(1.0, 1.0), Complex128(2.0, 2.0), Complex128(3.0, 3.0)],
            expect: [Complex128(0.0), Complex128(0.0), Complex128(0.0)]
        ),
        divTestStruct(
            dst: [Complex128.NaN, Complex128(1.0), Complex128.NaN, Complex128(1.0), Complex128(0.0), Complex128.NaN,
                Complex128(1.0), Complex128.NaN, Complex128(1.0), Complex128(0.0)],
            src: [Complex128(1.0), Complex128(1.0), Complex128.NaN, Complex128(1.0), Complex128(1.0), Complex128(1.0),
                Complex128(1.0), Complex128.NaN, Complex128(1.0), Complex128(1.0)],
            expect: [Complex128.NaN, Complex128(1.0), Complex128.NaN, Complex128(1.0), Complex128(0.0), Complex128.NaN,
                Complex128(1.0), Complex128.NaN, Complex128(1.0), Complex128(0.0)]
        ),
        divTestStruct(
            dst: [Complex128(inf), Complex128(4.0), Complex128.NaN, Complex128(inf), Complex128(9.0), Complex128(inf),
                Complex128(4.0), Complex128.NaN, Complex128(inf), Complex128(9.0)],
            src: [Complex128(inf), Complex128(4.0), Complex128.NaN, Complex128(inf), Complex128(3.0), Complex128(inf),
                Complex128(4.0), Complex128.NaN, Complex128(inf), Complex128(3.0)],
            expect: [Complex128.NaN, Complex128(1.0), Complex128.NaN, Complex128.NaN, Complex128(3.0), Complex128.NaN,
                Complex128(1.0), Complex128.NaN, Complex128.NaN, Complex128(3.0)]
        )
    ]
    for (j in 0..tests.size) {
        let sg_ln = 4 + j % 2
        let dg_ln = 4 + j % 3
        let src_ = guardVector(tests[j].src, src_gd, sg_ln)
        let dst_ = guardVector(tests[j].dst, dst_gd, dg_ln)
        let src = src_[sg_ln..src_.size - sg_ln]
        let dst = dst_[dg_ln..dst_.size - dg_ln]
        div(dst, src)
        for (i in 0..tests[j].expect.size) {
            @Expect(sameCmplxApprox(dst[i], tests[j].expect[i],tol))
        }
        @Expect(isValidGuard(src_, src_gd, sg_ln))
        @Expect(isValidGuard(dst_, dst_gd, dg_ln))
    }
}

class divToTestStruct {
    divToTestStruct(
        let dst!: Array<Complex128>,
        let x!: Array<Complex128>,
        let y!: Array<Complex128>,
        let expect!: Array<Complex128>
    ) {}
}

@Test
public func testDiTo() {
    const tol = 1e-15
    var dst_gd = Complex128(-1.0)
    var x_gd = Complex128(0.5)
    var y_gd = Complex128(0.25)
    let tests = [
        divToTestStruct(
            dst: [Complex128(1.0, -1.0)],
            x: [Complex128(1.0, 1.0)],
            y: [Complex128(1.0, 1.0)],
            expect: [Complex128(1.0)]
        ),
        divToTestStruct(
            dst: [Complex128(1.0, 0.0)],
            x: [Complex128.NaN],
            y: [Complex128.NaN],
            expect: [Complex128.NaN]
        ),
        divToTestStruct(
            dst: [Complex128(-2.0), Complex128(-2.0), Complex128(-2.0)],
            x: [Complex128(1.0, 1.0), Complex128(2.0, 1.0), Complex128(3.0, 1.0)],
            y: [Complex128(1.0, 1.0), Complex128(2.0, 1.0), Complex128(3.0, 1.0)],
            expect: [Complex128(1.0, 0.0), Complex128(1.0, 0.0), Complex128(1.0, 0.0)]
        ),
        divToTestStruct(
            dst: [Complex128(0.0), Complex128(0.0), Complex128(0.0)],
            x: [Complex128(2.0, 2.0), Complex128(4.0, 4.0), Complex128(6.0, 6.0)],
            y: [Complex128(1.0, 1.0), Complex128(2.0, 2.0), Complex128(3.0, 3.0)],
            expect: [Complex128(2.0, 0.0), Complex128(2.0, 0.0), Complex128(2.0, 0.0)]
        ),
        divToTestStruct(
            dst: [Complex128(-1.0), Complex128(-1.0), Complex128(-1.0)],
            x: [Complex128(0.0), Complex128(0.0), Complex128(0.0)],
            y: [Complex128(1.0, 1.0), Complex128(2.0, 1.0), Complex128(3.0, 1.0)],
            expect: [Complex128(0.0), Complex128(0.0), Complex128(0.0)]
        ),
        divToTestStruct(
            dst: [Complex128(inf), Complex128(inf), Complex128(inf), Complex128(inf), Complex128(inf), Complex128(inf),
                Complex128(inf), Complex128(inf), Complex128(inf), Complex128(inf)],
            x: [Complex128.NaN, Complex128(1.0), Complex128.NaN, Complex128(1.0), Complex128(0.0), Complex128.NaN,
                Complex128(1.0), Complex128.NaN, Complex128(1.0), Complex128(0.0)],
            y: [Complex128(1.0), Complex128(1.0), Complex128.NaN, Complex128(1.0), Complex128(1.0), Complex128(1.0),
                Complex128(1.0), Complex128.NaN, Complex128(1.0), Complex128(1.0)],
            expect: [Complex128.NaN, Complex128(1.0), Complex128.NaN, Complex128(1.0), Complex128(0.0), Complex128.NaN,
                Complex128(1.0), Complex128.NaN, Complex128(1.0), Complex128(0.0)]
        ),
        divToTestStruct(
            dst: [Complex128(0.0), Complex128(0.0), Complex128(0.0), Complex128(0.0), Complex128(0.0), Complex128(0.0),
                Complex128(0.0), Complex128(0.0), Complex128(0.0), Complex128(0.0)],
            x: [Complex128(inf), Complex128(4.0), Complex128.NaN, Complex128(inf), Complex128(9.0), Complex128(inf),
                Complex128(4.0), Complex128.NaN, Complex128(inf), Complex128(9.0)],
            y: [Complex128(inf), Complex128(4.0), Complex128.NaN, Complex128(inf), Complex128(3.0), Complex128(inf),
                Complex128(4.0), Complex128.NaN, Complex128(inf), Complex128(3.0)],
            expect: [Complex128.NaN, Complex128(1.0), Complex128.NaN, Complex128.NaN, Complex128(3.0), Complex128.NaN,
                Complex128(1.0), Complex128.NaN, Complex128.NaN, Complex128(3.0)]
        )
    ]
    for (j in 0..tests.size) {
        let xg_ln = 4 + j % 2
        let yg_ln = 4 + j % 3
        let y_ = guardVector(tests[j].y, y_gd, yg_ln)
        let x_ = guardVector(tests[j].x, x_gd, xg_ln)
        let y = y_[yg_ln..y_.size - yg_ln]
        let x = x_[xg_ln..x_.size - xg_ln]
        let dst_ = guardVector(tests[j].dst, dst_gd, xg_ln)
        let dst = dst_[xg_ln..dst_.size - xg_ln]
        let ret = divTo(dst, x, y)
        for (i in 0..tests[j].expect.size) {
            @Expect(sameCmplxApprox(ret[i], tests[j].expect[i],tol))
            @Expect(same(ret[i], dst[i]))
        }
        @Expect(isValidGuard(y_, y_gd, yg_ln))
        @Expect(isValidGuard(x_, x_gd, xg_ln))
        @Expect(isValidGuard(dst_, dst_gd, xg_ln))
    }
}

class sumTestStruct {
    sumTestStruct(
        let src!: Array<Complex128>,
        let expect!: Complex128
    ) {}
}

@Test
public func testSum() {
    var srcGd: Complex128 = Complex128(-0.1)
    let tests = [
        sumTestStruct(
            src: [],
            expect: Complex128(0.0)
        ),
        sumTestStruct(
            src: [Complex128(1.0)],
            expect: Complex128(1.0)
        ),
        sumTestStruct(
            src: [Complex128.NaN],
            expect: Complex128.NaN
        ),
        sumTestStruct(
            src: [Complex128(1.0, 1.0), Complex128(2.0, 2.0), Complex128(3.0, 3.0)],
            expect: Complex128(6.0, 6.0)
        ),
        sumTestStruct(
            src: [Complex128(1.0, 1.0), Complex128(-4.0), Complex128(3.0, -1.0)],
            expect: Complex128(0.0)
        ),
        sumTestStruct(
            src: [Complex128(1.0, -1.0), Complex128(2.0, 2.0), Complex128(3.0, -3.0), Complex128(4.0, 4.0)],
            expect: Complex128(10.0, 2.0)
        ),
        sumTestStruct(
            src: [Complex128(1.0), Complex128(1.0), Complex128.NaN, Complex128(1.0), Complex128(1.0)],
            expect: Complex128.NaN
        ),
        sumTestStruct(
            src: [Complex128(inf), Complex128(4.0), Complex128.NaN, Complex128(inf), Complex128(9.0)],
            expect: Complex128.NaN
        ),
        sumTestStruct(
            src: [Complex128(1.0, 1.0), Complex128(1.0, 1.0), Complex128(1.0, 1.0), Complex128(1.0, 1.0),
                Complex128(9.0, 9.0), Complex128(1.0, 1.0), Complex128(1.0, 1.0), Complex128(1.0, 1.0),
                Complex128(2.0, 2.0), Complex128(1.0, 1.0), Complex128(1.0, 1.0), Complex128(1.0, 1.0),
                Complex128(1.0, 1.0), Complex128(1.0, 1.0), Complex128(5.0, 5.0), Complex128(1.0, 1.0)],
            expect: Complex128(29.0, 29.0)
        ),
        sumTestStruct(
            src: [Complex128(1.0, 1.0), Complex128(1.0, 1.0), Complex128(1.0, 1.0), Complex128(1.0, 1.0),
                Complex128(9.0, 9.0), Complex128(1.0, 1.0), Complex128(1.0, 1.0), Complex128(1.0, 1.0),
                Complex128(2.0, 2.0), Complex128(1.0, 1.0), Complex128(1.0, 1.0), Complex128(1.0, 1.0),
                Complex128(1.0, 1.0), Complex128(1.0, 1.0), Complex128(5.0, 5.0), Complex128(11.0, 11.0),
                Complex128(1.0, 1.0), Complex128(1.0, 1.0), Complex128(1.0, 1.0), Complex128(9.0, 9.0),
                Complex128(1.0, 1.0), Complex128(1.0, 1.0), Complex128(1.0, 1.0), Complex128(2.0, 2.0),
                Complex128(1.0, 1.0), Complex128(1.0, 1.0), Complex128(1.0, 1.0), Complex128(1.0, 1.0),
                Complex128(1.0, 1.0), Complex128(5.0, 5.0), Complex128(1.0, 1.0)],
            expect: Complex128(67.0, 67.0)
        )
    ]
    for (j in 0..tests.size) {
        let gdLn = 4 + j % 2
        let gsrc = guardVector(tests[j].src, srcGd, gdLn)
        let src = gsrc[gdLn..gsrc.size - gdLn]
        let ret = sum(src)
        @Expect(same(ret, tests[j].expect))
        @Expect(isValidGuard(gsrc, srcGd, gdLn))
    }
}
