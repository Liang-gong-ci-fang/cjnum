package cjnum.internal.asm.c128

import cjnum.cmplxs.cscalar.*
import cjnum.complex.*
import std.unittest.*
import std.unittest.testmacro.*

class axpyTestStruct {
    axpyTestStruct(
        let incX!: Int64,
        let incY!: Int64,
        let incDst!: Int64,
        let ix!: Int64,
        let iy!: Int64,
        let idst!: Int64,
        let a!: Complex128,
        let dst!: Array<Complex128>,
        let x!: Array<Complex128>,
        let y!: Array<Complex128>,
        let ex!: Array<Complex128>
    ) {}
}

let axpyTests = [
    axpyTestStruct(incX: 2, incY: 2, incDst: 3, ix: 0, iy: 0, idst: 0, a: Complex128(1.0, 1.0),
        dst: [Complex128(5.0, 0.0)], x: [Complex128(1.0, 0.0)], y: [Complex128(0.0, 1.0)], ex: [Complex128(1.0, 2.0)]),
    axpyTestStruct(incX: 2, incY: 2, incDst: 3, ix: 0, iy: 0, idst: 0, a: Complex128(1.0, 2.0),
        dst: [Complex128(0.0, 0.0), Complex128(0.0, 0.0), Complex128(0.0, 0.0)],
        x: [Complex128(0.0, 0.0), Complex128(0.0, 0.0), Complex128(0.0, 0.0)],
        y: [Complex128(1.0, 0.0), Complex128(1.0, 0.0), Complex128(1.0, 0.0)],
        ex: [Complex128(1.0, 0.0), Complex128(1.0, 0.0), Complex128(1.0, 0.0)]),
    axpyTestStruct(incX: 2, incY: 2, incDst: 3, ix: 0, iy: 0, idst: 0, a: Complex128(1.0, 2.0),
        dst: [Complex128(0.0, 0.0), Complex128(0.0, 0.0), Complex128(0.0, 0.0)],
        x: [Complex128(0.0, 0.0), Complex128(0.0, 0.0)],
        y: [Complex128(1.0, 0.0), Complex128(1.0, 0.0), Complex128(1.0, 0.0)],
        ex: [Complex128(1.0, 0.0), Complex128(1.0, 0.0)]),
    axpyTestStruct(incX: 2, incY: 2, incDst: 3, ix: 0, iy: 0, idst: 0, a: Complex128(1.0, 2.0),
        dst: [Complex128(0.0, 1.0), Complex128(0.0, 1.0), Complex128(0.0, 1.0)],
        x: [Complex128(0.0, 1.0), Complex128(0.0, 1.0), Complex128(0.0, 1.0)],
        y: [Complex128(1.0, 0.0), Complex128(2.0, 0.0), Complex128(1.0, 0.0)],
        ex: [Complex128(-1.0, 1.0), Complex128(0.0, 1.0), Complex128(-1.0, 1.0)]),
    axpyTestStruct(incX: 2, incY: 2, incDst: 3, ix: 0, iy: 0, idst: 0, a: Complex128(0.0, -1.0),
        dst: [Complex128(0.0, 1.0), Complex128(0.0, 1.0), Complex128(0.0, 1.0)],
        x: [Complex128(0.0, 1.0), Complex128(0.0, 1.0), Complex128(0.0, 1.0)],
        y: [Complex128(1.0, 0.0), Complex128(2.0, 0.0), Complex128(1.0, 0.0)],
        ex: [Complex128(2.0, 0.0), Complex128(3.0, 0.0), Complex128(2.0, 0.0)]),
    axpyTestStruct(incX: 2, incY: 2, incDst: 3, ix: 0, iy: 0, idst: 0, a: Complex128(0.0, -1.0),
        dst: [Complex128(0.0, 1.0), Complex128(0.0, 1.0), Complex128(0.0, 1.0)],
        x: [Complex128(0.0, 1.0), Complex128(0.0, 1.0), Complex128(0.0, 1.0)],
        y: [Complex128(1.0, 0.0), Complex128(2.0, 0.0), Complex128(1.0, 0.0)],
        ex: [Complex128(2.0, 0.0), Complex128(3.0, 0.0), Complex128(2.0, 0.0)]),
    axpyTestStruct(incX: 2, incY: 4, incDst: 3, ix: 0, iy: 0, idst: 0, a: Complex128(-2.0, 0.0),
        dst: [Complex128(0.0, 1.0), Complex128(0.0, 1.0), Complex128(0.0, 1.0), Complex128(0.0, 1.0),
            Complex128(0.0, 1.0)],
        x: [Complex128(2.0, 1.0), Complex128(2.0, 1.0), Complex128(2.0, 1.0), Complex128(2.0, 1.0), Complex128(2.0, 1.0)],
        y: [Complex128(1.0, 0.0), Complex128(1.0, 0.0), Complex128(2.0, 0.0), Complex128(1.0, 0.0), Complex128(1.0, 0.0)],
        ex: [Complex128(-3.0, -2.0), Complex128(-3.0, -2.0), Complex128(-2.0, -2.0), Complex128(-3.0, -2.0),
            Complex128(-3.0, -2.0)]),
    axpyTestStruct(incX: 2, incY: 2, incDst: 3, ix: 0, iy: 0, idst: 0, a: Complex128(1.0, -1.0),
        dst: Array<Complex128>(10, repeat: Complex128(0.0, 0.0)), x: Array<Complex128>(10, repeat: Complex128(0.0, 1.0)),
        y: [Complex128(1.0, 0.0), Complex128(1.0, 0.0), Complex128(2.0, 0.0), Complex128(1.0, 0.0), Complex128(1.0, 0.0),
            Complex128(1.0, 0.0), Complex128(1.0, 0.0), Complex128(2.0, 0.0), Complex128(1.0, 0.0), Complex128(1.0, 0.0)],
        ex: [Complex128(2.0, 1.0), Complex128(2.0, 1.0), Complex128(3.0, 1.0), Complex128(2.0, 1.0), Complex128(2.0, 1.0
        ), Complex128(2.0, 1.0), Complex128(2.0, 1.0), Complex128(3.0, 1.0), Complex128(2.0, 1.0), Complex128(2.0, 1.0)]
    ),
    axpyTestStruct(incX: 2, incY: 2, incDst: 3, ix: 0, iy: 0, idst: 0, a: Complex128(1.0, -1.0),
        dst: Array<Complex128>(10, repeat: Complex128(0.0, 0.0)), x: Array<Complex128>(10, repeat: Complex128(0.0, 1.0)),
        y: [Complex128(1.0, 0.0), Complex128(1.0, 0.0), Complex128(2.0, 0.0), Complex128(1.0, 0.0), Complex128(1.0, 0.0),
            Complex128(1.0, 0.0), Complex128(1.0, 0.0), Complex128(2.0, 0.0), Complex128(1.0, 0.0), Complex128(1.0, 0.0)],
        ex: [Complex128(2.0, 1.0), Complex128(2.0, 1.0), Complex128(3.0, 1.0), Complex128(2.0, 1.0), Complex128(2.0, 1.0
        ), Complex128(2.0, 1.0), Complex128(2.0, 1.0), Complex128(3.0, 1.0), Complex128(2.0, 1.0), Complex128(2.0, 1.0)]
    ),
    axpyTestStruct(incX: -2, incY: -2, incDst: -3, ix: 18, iy: 18, idst: 27, a: Complex128(1.0, -1.0),
        dst: Array<Complex128>(10, repeat: Complex128(0.0, 0.0)), x: Array<Complex128>(10, repeat: Complex128(0.0, 1.0)),
        y: [Complex128(1.0, 0.0), Complex128(1.0, 0.0), Complex128(2.0, 0.0), Complex128(1.0, 0.0), Complex128(1.0, 0.0),
            Complex128(1.0, 0.0), Complex128(1.0, 0.0), Complex128(2.0, 0.0), Complex128(1.0, 0.0), Complex128(1.0, 0.0)],
        ex: [Complex128(2.0, 1.0), Complex128(2.0, 1.0), Complex128(3.0, 1.0), Complex128(2.0, 1.0), Complex128(2.0, 1.0
        ), Complex128(2.0, 1.0), Complex128(2.0, 1.0), Complex128(3.0, 1.0), Complex128(2.0, 1.0), Complex128(2.0, 1.0)]
    ),
    axpyTestStruct(incX: -2, incY: 2, incDst: -3, ix: 18, iy: 0, idst: 27, a: Complex128(1.0, -1.0),
        dst: Array<Complex128>(10, repeat: Complex128(0.0, 0.0)), x: Array<Complex128>(10, repeat: Complex128(0.0, 1.0)),
        y: [Complex128(1.0, 0.0), Complex128(1.0, 0.0), Complex128(2.0, 0.0), Complex128(1.0, 0.0), Complex128(1.0, 0.0),
            Complex128(1.0, 0.0), Complex128(1.0, 0.0), Complex128(2.0, 0.0), Complex128(1.0, 0.0), Complex128(1.0, 0.0)],
        ex: [Complex128(2.0, 1.0), Complex128(2.0, 1.0), Complex128(3.0, 1.0), Complex128(2.0, 1.0), Complex128(2.0, 1.0
        ), Complex128(2.0, 1.0), Complex128(2.0, 1.0), Complex128(3.0, 1.0), Complex128(2.0, 1.0), Complex128(2.0, 1.0)]
    )
]

@Test
public func testAxpyUnitary() {
    let xGdVal: Complex128 = Complex128(1.0, 0.0)
    let yGdVal: Complex128 = Complex128(1.0, 0.0)
    for (cas in 0..axpyTests.size) {
        let xgLn = 4 + cas % 2
        let ygLn = 4 + cas % 3
        let axpy_x = guardVector(axpyTests[cas].x, xGdVal, xgLn)
        let axpy_y = guardVector(axpyTests[cas].y, yGdVal, ygLn)
        let x = axpy_x[xgLn..(axpy_x.size) - xgLn]
        let y = axpy_y[ygLn..(axpy_y.size) - ygLn]
        axpyUnitary(axpyTests[cas].a, x, y)
        for (i in 0..axpyTests[cas].ex.size) {
            @Expect(y[i]==axpyTests[cas].ex[i])
        }
        @Expect(isValidGuard(axpy_x, xGdVal, xgLn))
        @Expect(isValidGuard(axpy_y, yGdVal, ygLn))
    }
}

@Test
public func testAxpyUnitaryTo() {
    let xGdVal: Complex128 = Complex128(1.0, 0.0)
    let yGdVal: Complex128 = Complex128(1.0, 0.0)
    let dstGdVal: Complex128 = Complex128(0.0, 0.0)
    for (cas in 0..axpyTests.size) {
        let xgLn = 4 + cas % 2
        let ygLn = 4 + cas % 3
        let axpy_x = guardVector(axpyTests[cas].x, xGdVal, xgLn)
        let axpy_y = guardVector(axpyTests[cas].y, yGdVal, ygLn)
        let axpy_dst = guardVector(axpyTests[cas].dst, dstGdVal, xgLn)
        let x = axpy_x[xgLn..(axpy_x.size) - xgLn]
        let y = axpy_y[ygLn..(axpy_y.size) - ygLn]
        let dst = axpy_dst[xgLn..(axpy_dst.size) - xgLn]
        axpyUnitaryTo(dst, axpyTests[cas].a, x, y)
        for (i in 0..axpyTests[cas].ex.size) {
            @Expect(dst[i]==axpyTests[cas].ex[i])
        }
        @Expect(isValidGuard(axpy_x, xGdVal, xgLn))
        @Expect(isValidGuard(axpy_y, yGdVal, ygLn))
        @Expect(isValidGuard(axpy_dst, dstGdVal, xgLn))
    }
}

@Test
public func testAxpyInc() {
    let xGdVal: Complex128 = Complex128(1.0)
    let yGdVal: Complex128 = Complex128(1.0)
    for (cas in 0..axpyTests.size) {
        let xgLn = 4 + cas % 2
        let ygLn = 4 + cas % 3
        let axpy_x = guardIncVector(axpyTests[cas].x, xGdVal, axpyTests[cas].incX, xgLn)
        let axpy_y = guardIncVector(axpyTests[cas].y, yGdVal, axpyTests[cas].incY, ygLn)
        let x = axpy_x[xgLn..(axpy_x.size - xgLn)]
        let y = axpy_y[ygLn..(axpy_y.size - ygLn)]
        axpyInc(axpyTests[cas].a, x, y, axpyTests[cas].ex.size, axpyTests[cas].incX, axpyTests[cas].incY,
            axpyTests[cas].ix, axpyTests[cas].iy)
        for (i in 0..axpyTests[cas].ex.size) {
            @Expect(y[axpyTests[cas].iy+i*axpyTests[cas].incY]==axpyTests[cas].ex[i])
        }
        @Expect(checkValidIncGuard(axpy_x, xGdVal,axpyTests[cas].incX, xgLn))
        @Expect(checkValidIncGuard(axpy_y, yGdVal,axpyTests[cas].incY, ygLn))
    }
}

@Test
public func testAxpyIncTo() {
    let xGdVal: Complex128 = Complex128(1.0, 0.0)
    let yGdVal: Complex128 = Complex128(1.0, 0.0)
    let dstGdVal: Complex128 = Complex128(0.0, 0.0)
    for (cas in 0..axpyTests.size) {
        let xgLn = 4 + cas % 2
        let ygLn = 4 + cas % 3
        let axpy_x = guardIncVector(axpyTests[cas].x, xGdVal, axpyTests[cas].incX, xgLn)
        let axpy_y = guardIncVector(axpyTests[cas].y, yGdVal, axpyTests[cas].incY, ygLn)
        let axpy_dst = guardIncVector(axpyTests[cas].dst, dstGdVal, axpyTests[cas].incDst, xgLn)
        let x = axpy_x[xgLn..(axpy_x.size) - xgLn]
        let y = axpy_y[ygLn..(axpy_y.size) - ygLn]
        let dst = axpy_dst[xgLn..(axpy_dst.size) - xgLn]
        axpyIncTo(dst, axpyTests[cas].incDst, axpyTests[cas].idst, axpyTests[cas].a, x, y, axpyTests[cas].ex.size,
            axpyTests[cas].incX, axpyTests[cas].incY, axpyTests[cas].ix, axpyTests[cas].iy)
        for (i in 0..axpyTests[cas].ex.size) {
            @Expect(dst[axpyTests[cas].idst+i*axpyTests[cas].incDst]==axpyTests[cas].ex[i])
        }
        @Expect(checkValidIncGuard(axpy_x, xGdVal,axpyTests[cas].incX, xgLn))
        @Expect(checkValidIncGuard(axpy_y, yGdVal, axpyTests[cas].incY,ygLn))
        @Expect(checkValidIncGuard(axpy_dst, dstGdVal,axpyTests[cas].incDst, xgLn))
    }
}

class dscalTestStruct {
    dscalTestStruct(
        let alpha!: Float64,
        let x!: Array<Complex128>,
        let want!: Array<Complex128>
    ) {}
}

let dscalTests = [
    dscalTestStruct(alpha: 0.0, x: [Complex128(0.0, 0.0)], want: [Complex128(0.0, 0.0)]),
    dscalTestStruct(alpha: 1.0, x: [Complex128(1.0, 2.0)], want: [Complex128(1.0, 2.0)]),
    dscalTestStruct(alpha: 2.0, x: [Complex128(1.0, 2.0)], want: [Complex128(2.0, 4.0)]),
    dscalTestStruct(
        alpha: 2.0,
        x: [Complex128(1.0, 2.0), Complex128(3.0, 5.0), Complex128(6.0, 11.0), Complex128(12.0, -23.0)],
        want: [Complex128(2.0, 4.0), Complex128(6.0, 10.0), Complex128(12.0, 22.0), Complex128(24.0, -46.0)]
    ),
    dscalTestStruct(
        alpha: 3.0,
        x: [Complex128(1.0, 2.0), Complex128(5.0, 4.0), Complex128(3.0, 6.0), Complex128(8.0, 12.0),
            Complex128(-3.0, -2.0), Complex128(-5.0, 5.0)],
        want: [Complex128(3.0, 6.0), Complex128(15.0, 12.0), Complex128(9.0, 18.0), Complex128(24.0, 36.0),
            Complex128(-9.0, -6.0), Complex128(-15.0, 15.0)]
    ),
    dscalTestStruct(
        alpha: 5.0,
        x: [Complex128(1.0, 2.0), Complex128(5.0, 4.0), Complex128(3.0, 6.0), Complex128(8.0, 12.0),
            Complex128(-3.0, -2.0), Complex128(-5.0, 5.0), Complex128(1.0, 2.0), Complex128(5.0, 4.0),
            Complex128(3.0, 6.0), Complex128(8.0, 12.0), Complex128(-3.0, -2.0), Complex128(-5.0, 5.0)],
        want: [Complex128(5.0, 10.0), Complex128(25.0, 20.0), Complex128(15.0, 30.0), Complex128(40.0, 60.0),
            Complex128(-15.0, -10.0), Complex128(-25.0, 25.0), Complex128(5.0, 10.0), Complex128(25.0, 20.0),
            Complex128(15.0, 30.0), Complex128(40.0, 60.0), Complex128(-15.0, -10.0), Complex128(-25.0, 25.0)]
    )
]

@Test
public func testDscalUnitary() {
    let xGdVal: Complex128 = Complex128(-0.5, 0.0)
    for (test in dscalTests) {
        for (align in align1) {
            let xgLn = 4 + align
            let xg = guardVector(test.x, xGdVal, xgLn)
            let x = xg[xgLn..(xg.size) - xgLn]
            dscalUnitary(test.alpha, x)
            for (i in 0..test.want.size) {
                @Expect(same(x[i], test.want[i]))
            }
            @Expect(isValidGuard(xg, xGdVal, xgLn))
        }
    }
}

@Test
public func testDscalInc() {
    let xGdVal: Complex128 = Complex128(-0.5)
    let gdLn = 4
    for (test in dscalTests) {
        let n = test.x.size
        for (incX in [1, 2, 3, 4, 7, 10]) {
            let xg = guardIncVector(test.x, xGdVal, incX, gdLn)
            let x = xg[gdLn..(xg.size - gdLn)]
            dscalInc(Complex128(test.alpha), x, n, incX)
            for (i in 0..test.want.size) {
                @Expect(same(x[i*incX], test.want[i]))
            }
            @Expect(checkValidIncGuard(xg, xGdVal,incX, gdLn))
        }
    }
}

class scalTestStruct {
    scalTestStruct(
        let alpha!: Complex128,
        let x!: Array<Complex128>,
        let want!: Array<Complex128>
    ) {}
}

let scalTests = [
    scalTestStruct(alpha: Complex128(0.0, 0.0), x: [], want: []),
    scalTestStruct(alpha: Complex128(1.0, 1.0), x: [Complex128(1.0, 2.0)], want: [Complex128(-1.0, 3.0)]),
    scalTestStruct(alpha: Complex128(2.0, 3.0), x: [Complex128(1.0, 2.0)], want: [Complex128(-4.0, 7.0)]),
    scalTestStruct(alpha: Complex128(2.0, -4.0), x: [Complex128(1.0, 2.0)], want: [Complex128(10.0, 0.0)]),
    scalTestStruct(
        alpha: Complex128(2.0, 8.0),
        x: [Complex128(1.0, 2.0), Complex128(5.0, 4.0), Complex128(3.0, 6.0), Complex128(8.0, 12.0),
            Complex128(-3.0, -2.0), Complex128(-5.0, 5.0)],
        want: [Complex128(-14.0, 12.0), Complex128(-22.0, 48.0), Complex128(-42.0, 36.0), Complex128(-80.0, 88.0),
            Complex128(10.0, -28.0), Complex128(-50.0, -30.0)]
    ),
    scalTestStruct(
        alpha: Complex128(5.0, -10.0),
        x: [Complex128(1.0, 2.0), Complex128(5.0, 4.0), Complex128(3.0, 6.0), Complex128(8.0, 12.0),
            Complex128(-3.0, -2.0), Complex128(-5.0, 5.0), Complex128(1.0, 2.0), Complex128(5.0, 4.0),
            Complex128(3.0, 6.0), Complex128(8.0, 12.0), Complex128(-3.0, -2.0), Complex128(-5.0, 5.0)],
        want: [Complex128(25.0, 0.0), Complex128(65.0, -30.0), Complex128(75.0, 0.0), Complex128(160.0, -20.0),
            Complex128(-35.0, 20.0), Complex128(25.0, 75.0), Complex128(25.0, 0.0), Complex128(65.0, -30.0),
            Complex128(75.0, 0.0), Complex128(160.0, -20.0), Complex128(-35.0, 20.0), Complex128(25.0, 75.0)]
    )
]

@Test
public func testScalUnitary() {
    let xGdVal: Complex128 = Complex128(-0.5, 0.0)
    for (test in scalTests) {
        for (align in align1) {
            let xgLn = 4 + align
            let xg = guardVector(test.x, xGdVal, xgLn)
            let x = xg[xgLn..(xg.size) - xgLn]
            scalUnitary(test.alpha, x)
            for (i in 0..test.want.size) {
                @Expect(same(x[i], test.want[i]))
            }
            @Expect(isValidGuard(xg, xGdVal, xgLn))
        }
    }
}

@Test
public func testScalInc() {
    let xGdVal: Complex128 = Complex128(-0.5)
    let gdLn = 4
    for (test in scalTests) {
        let n = test.x.size
        for (incX in [1, 2, 3, 4, 7, 10]) {
            let xg = guardIncVector(test.x, xGdVal, incX, gdLn)
            let x = xg[gdLn..(xg.size - gdLn)]
            scalInc(test.alpha, x, n, incX)
            for (i in 0..test.want.size) {
                @Expect(same(x[i*incX], test.want[i]))
            }
            @Expect(checkValidIncGuard(xg, xGdVal,incX, gdLn))
        }
    }
}
