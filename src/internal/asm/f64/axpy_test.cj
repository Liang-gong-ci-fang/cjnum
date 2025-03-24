package cjnum.internal.asm.f64

import std.math.*
import std.unittest.*
import std.unittest.testmacro.*
import cjnum.floats.scalar.*

class axpyTestStruct {
    axpyTestStruct(
        let alpha: Float64,
        let x: Array<Float64>,
        let y: Array<Float64>,
        let want: Array<Float64>,
        let wantRev: Array<Float64>
    ) {}
}

let axpyTests = [
    axpyTestStruct(0.0, Array<Float64>(), Array<Float64>(), Array<Float64>(), Array<Float64>()),
    axpyTestStruct(0.0, Array<Float64>([2.0]), Array<Float64>([-3.0]), Array<Float64>([-3.0]), Array<Float64>([-3.0])),
    axpyTestStruct(1.0, Array<Float64>([2.0]), Array<Float64>([-3.0]), Array<Float64>([-1.0]), Array<Float64>([-1.0])),
    axpyTestStruct(3.0, Array<Float64>([2.0]), Array<Float64>([-3.0]), Array<Float64>([3.0]), Array<Float64>([3.0])),
    axpyTestStruct(-3.0, Array<Float64>([2.0]), Array<Float64>([-3.0]), Array<Float64>([-9.0]), Array<Float64>([-9.0])),
    axpyTestStruct(1.0, Array<Float64>([1.0, 5.0]), Array<Float64>([2.0, -3.0]), Array<Float64>([3.0, 2.0]),
        Array<Float64>([7.0, -2.0])),
    axpyTestStruct(1.0, Array<Float64>([2.0, 3.0, 4.0]), Array<Float64>([-3.0, -2.0, -1.0]),
        Array<Float64>([-1.0, 1.0, 3.0]), Array<Float64>([1.0, 1.0, 1.0])),
    axpyTestStruct(0.0, Array<Float64>([0.0, 0.0, 1.0, 1.0, 2.0, -3.0, -4.0]),
        Array<Float64>([0.0, 1.0, 0.0, 3.0, -4.0, 5.0, -6.0]), Array<Float64>([0.0, 1.0, 0.0, 3.0, -4.0, 5.0, -6.0]),
        Array<Float64>([0.0, 1.0, 0.0, 3.0, -4.0, 5.0, -6.0])),
    axpyTestStruct(1.0, Array<Float64>([0.0, 0.0, 1.0, 1.0, 2.0, -3.0, -4.0]),
        Array<Float64>([0.0, 1.0, 0.0, 3.0, -4.0, 5.0, -6.0]), Array<Float64>([0.0, 1.0, 1.0, 4.0, -2.0, 2.0, -10.0]),
        Array<Float64>([-4.0, -2.0, 2.0, 4.0, -3.0, 5.0, -6.0])),
    axpyTestStruct(3.0, Array<Float64>([0.0, 0.0, 1.0, 1.0, 2.0, -3.0, -4.0]),
        Array<Float64>([0.0, 1.0, 0.0, 3.0, -4.0, 5.0, -6.0]), Array<Float64>([0.0, 1.0, 3.0, 6.0, 2.0, -4.0, -18.0]),
        Array<Float64>([-12.0, -8.0, 6.0, 6.0, -1.0, 5.0, -6.0])),
    axpyTestStruct(-3.0, Array<Float64>([0.0, 0.0, 1.0, 1.0, 2.0, -3.0, -4.0, 0.0, 0.0, 1.0, 1.0, 2.0, -3.0, -4.0]),
        Array<Float64>([0.0, 1.0, 0.0, 3.0, -4.0, 5.0, -6.0, 0.0, 1.0, 0.0, 3.0, -4.0, 5.0, -6.0]),
        Array<Float64>([0.0, 1.0, -3.0, 0.0, -10.0, 14.0, 6.0, 0.0, 1.0, -3.0, 0.0, -10.0, 14.0, 6.0]),
        Array<Float64>([12.0, 10.0, -6.0, 0.0, -7.0, 5.0, -6.0, 12.0, 10.0, -6.0, 0.0, -7.0, 5.0, -6.0])),
    axpyTestStruct(-5.0, Array<Float64>([0.0, 0.0, 1.0, 1.0, 2.0, -3.0, -4.0, 5.0, 1.0, 2.0, -3.0, -4.0, 5.0]),
        Array<Float64>([0.0, 1.0, 0.0, 3.0, -4.0, 5.0, -6.0, 7.0, 3.0, -4.0, 5.0, -6.0, 7.0]),
        Array<Float64>([0.0, 1.0, -5.0, -2.0, -14.0, 20.0, 14.0, -18.0, -2.0, -14.0, 20.0, 14.0, -18.0]),
        Array<Float64>([-25.0, 21.0, 15.0, -7.0, -9.0, -20.0, 14.0, 22.0, -7.0, -9.0, 0.0, -6.0, 7.0]))
]

@Test
public func testAxpyUnitary() {
    const xGdVal: Float64 = -1.0
    const yGdVal: Float64 = 0.5
    for (i in 0..axpyTests.size) {
        for (align in align2) {
            let xgLn = 4 + align.x
            let ygLn = 4 + align.y
            let xg = guardVector(axpyTests[i].x, xGdVal, xgLn)
            let yg = guardVector(axpyTests[i].y, yGdVal, ygLn)
            let x = xg[xgLn..(xg.size - xgLn)]
            let y = yg[ygLn..(yg.size - ygLn)]

            axpyUnitary(axpyTests[i].alpha, x, y)
            for (j in 0..axpyTests[i].want.size) {
                @Expect(same(y[j], axpyTests[i].want[j]))
            }
            @Expect(isValidGuard(xg, xGdVal, xgLn))
            @Expect(isValidGuard(yg, yGdVal, ygLn))
            @Expect(equalStrided(axpyTests[i].x, x, 1))
        }
    }
}

@Test
public func testAxpyUnitaryTo() {
    const dstGdVal: Float64 = 1.0
    const xGdVal: Float64 = -1.0
    const yGdVal: Float64 = 0.5
    for (i in 0..axpyTests.size) {
        for (align in align3) {
            let dgLn = 4 + align.dst
            let xgLn = 4 + align.x
            let ygLn = 4 + align.y
            let dstOrig = Array<Float64>(axpyTests[i].x.size, item: 0.0)
            let xg = guardVector(axpyTests[i].x, xGdVal, xgLn)
            let yg = guardVector(axpyTests[i].y, yGdVal, ygLn)
            let dstg = guardVector(dstOrig, dstGdVal, dgLn)
            let x = xg[xgLn..(xg.size - xgLn)]
            let y = yg[ygLn..(yg.size - ygLn)]
            let dst = dstg[dgLn..(dstg.size - dgLn)]

            axpyUnitaryTo(dst, axpyTests[i].alpha, x, y)
            for (j in 0..axpyTests[i].want.size) {
                @Expect(same(dst[j], axpyTests[i].want[j]))
            }
            @Expect(isValidGuard(xg, xGdVal, xgLn))
            @Expect(isValidGuard(yg, yGdVal, ygLn))
            @Expect(isValidGuard(dstg, dstGdVal, dgLn))
            @Expect(equalStrided(axpyTests[i].x, x, 1))
            @Expect(equalStrided(axpyTests[i].y, y, 1))
        }
    }
}

@Test
public func testAxpyInc() {
    const xGdVal: Float64 = -1.0
    const yGdVal: Float64 = 0.5
    let gdLn = 4
    for (i in 0..axpyTests.size) {
        let n = axpyTests[i].x.size
        for (inc in newIncSet([-7, -4, -3, -2, -1, 1, 2, 3, 4, 7])) {
            var ix: Int64 = 0
            var iy: Int64 = 0
            if (inc.x < 0) {
                ix = (-n + 1) * inc.x
            }
            if (inc.y < 0) {
                iy = (-n + 1) * inc.y
            }
            let xg = guardIncVector(axpyTests[i].x, xGdVal, inc.x, gdLn)
            let yg = guardIncVector(axpyTests[i].y, yGdVal, inc.y, gdLn)
            let x = xg[gdLn..(xg.size - gdLn)]
            let y = yg[gdLn..(yg.size - gdLn)]

            axpyInc(axpyTests[i].alpha, x, y, n, inc.x, inc.y, ix, iy)

            var want = axpyTests[i].want
            if (inc.x * inc.y < 0) {
                want = axpyTests[i].wantRev
            }
            inc.y = abs(inc.y)
            for (j in 0..want.size) {
                @Expect(same(y[j * inc.y], want[j]))
            }
            @Expect(equalStrided(axpyTests[i].x, x, inc.x))
            @Expect(checkValidIncGuard(xg, xGdVal, inc.x, gdLn))
            @Expect(checkValidIncGuard(yg, yGdVal, inc.y, gdLn))
        }
    }
}

@Test
public func testAxpyIncTo() {
    const dstGdVal: Float64 = 1.0
    const xGdVal: Float64 = -1.0
    const yGdVal: Float64 = 0.5
    var want = Array<Float64>()
    let gdLn = 4
    for (i in 0..axpyTests.size) {
        let n = axpyTests[i].x.size
        for (inc in newIncToSet([-7, -4, -3, -2, -1, 1, 2, 3, 4, 7])) {
            var ix: Int64 = 0
            var iy: Int64 = 0
            var idst: Int64 = 0
            if (inc.x < 0) {
                ix = (-n + 1) * inc.x
            }
            if (inc.y < 0) {
                iy = (-n + 1) * inc.y
            }
            if (inc.dst < 0) {
                idst = (-n + 1) * inc.dst
            }

            let dstOrig = Array<Float64>(axpyTests[i].want.size, item: 0.0)
            let xg = guardIncVector(axpyTests[i].x, xGdVal, inc.x, gdLn)
            let yg = guardIncVector(axpyTests[i].y, yGdVal, inc.y, gdLn)
            let dstg = guardIncVector(dstOrig, dstGdVal, inc.dst, gdLn)
            let x = xg[gdLn..(xg.size - gdLn)]
            let y = yg[gdLn..(yg.size - gdLn)]
            let dst = dstg[gdLn..(dstg.size - gdLn)]

            axpyIncTo(dst, inc.dst, idst, axpyTests[i].alpha, x, y, n, inc.x, inc.y, ix, iy)
            want = axpyTests[i].want
            if (inc.x * inc.y < 0) {
                want = axpyTests[i].wantRev
            }
            var iW = 0
            var incW = 1
            if (inc.y * inc.dst < 0) {
                iW = want.size - 1
                incW = -1
            }
            inc.dst = abs(inc.dst)
            for (j in 0..want.size) {
                @Expect(same(dst[j * inc.dst], want[iW + j * incW]))
            }

            @Expect(checkValidIncGuard(xg, xGdVal, inc.x, gdLn))
            @Expect(checkValidIncGuard(yg, yGdVal, inc.y, gdLn))
            @Expect(checkValidIncGuard(dstg, dstGdVal, inc.dst, gdLn))
            @Expect(equalStrided(axpyTests[i].x, x, inc.x))
            @Expect(equalStrided(axpyTests[i].y, y, inc.y))
        }
    }
}
