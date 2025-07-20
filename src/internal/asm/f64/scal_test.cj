package cjnum.internal.asm.f64

import std.unittest.*
import std.unittest.testmacro.*
import cjnum.floats.scalar.*
import std.random.*

class scalTestStruct {
    scalTestStruct(
        let alpha: Float64,
        let x: Array<Float64>,
        let want: Array<Float64>
    ) {}
}

let scalTests = [
    scalTestStruct(0.0, Array<Float64>(), Array<Float64>()),
    scalTestStruct(0.0, [1.0], [0.0]),
    scalTestStruct(1.0, [1.0], [1.0]),
    scalTestStruct(2.0, [1.0, -2.0], [2.0, -4.0]),
    scalTestStruct(2.0, [1.0, -2.0, 3.0], [2.0, -4.0, 6.0]),
    scalTestStruct(2.0, [1.0, -2.0, 3.0, 4.0], [2.0, -4.0, 6.0, 8.0]),
    scalTestStruct(2.0, [1.0, -2.0, 3.0, 4.0, -5.0], [2.0, -4.0, 6.0, 8.0, -10.0]),
    scalTestStruct(2.0, [0.0, 1.0, -2.0, 3.0, 4.0, -5.0, 6.0, -7.0],
        [0.0, 2.0, -4.0, 6.0, 8.0, -10.0, 12.0, -14.0]),
    scalTestStruct(2.0, [0.0, 1.0, -2.0, 3.0, 4.0, -5.0, 6.0, -7.0, 8.0],
        [0.0, 2.0, -4.0, 6.0, 8.0, -10.0, 12.0, -14.0, 16.0]),
    scalTestStruct(2.0, [0.0, 1.0, -2.0, 3.0, 4.0, -5.0, 6.0, -7.0, 8.0, 9.0],
        [0.0, 2.0, -4.0, 6.0, 8.0, -10.0, 12.0, -14.0, 16.0, 18.0]),
    scalTestStruct(3.0, [0.0, 1.0, -2.0, 3.0, 4.0, -5.0, 6.0, -7.0, 8.0, 9.0, 12.0],
        [0.0, 3.0, -6.0, 9.0, 12.0, -15.0, 18.0, -21.0, 24.0, 27.0, 36.0])
]

@Test
public func testScalUnitary() {
    const xGdVal = -0.5

    for (i in 0..scalTests.size) {
        let test = scalTests[i]
        for (align in align1) {
            let xgLn: Int64 = 4 + align
            let xg = guardVector(test.x, xGdVal, xgLn)
            let x = xg[xgLn..xg.size - xgLn]

            scalUnitary(test.alpha, x)

            for (i in 0..test.want.size) {
                @Expect(same(x[i], test.want[i]))
            }
            @Expect(isValidGuard(xg, xGdVal, xgLn))
        }
    }
}

@Test
public func testScalUnitaryTo() {
    const xGdVal = -1.0
    const dstGdVal = 0.5
    let rng = Random(42)
    for (i in 0..scalTests.size) {
        let test = scalTests[i]
        let n = test.x.size
        for (align in align2) {
            let xgLn: Int64 = 4 + align.x
            let dgLn: Int64 = 4 + align.y
            let xg = guardVector(test.x, xGdVal, xgLn)
            let dg = guardVector(randSlice(n, 1, rng), dstGdVal, dgLn)
            let x = xg[xgLn..xg.size - xgLn]
            let dst = dg[dgLn..dg.size - dgLn]

            scalUnitaryTo(dst, test.alpha, x)

            for (i in 0..test.want.size) {
                @Expect(same(dst[i], test.want[i]))
            }
            @Expect(isValidGuard(xg, xGdVal, xgLn))
            @Expect(isValidGuard(dg, dstGdVal, dgLn))
            @Expect(equalStrided(test.x, x, 1))
        }
    }
}

@Test
public func testScalInc() {
    const xGdVal = -0.5
    let gdln = 4
    for (i in 0..scalTests.size) {
        let test = scalTests[i]
        let n = test.x.size
        for (incX in [1, 2, 3, 4, 7, 10]) {
            let xg = guardIncVector(test.x, xGdVal, incX, gdln)
            let x = xg[gdln..xg.size - gdln]

            scalInc(test.alpha, x, n, incX)

            for (i in 0..test.want.size) {
                @Expect(same(x[i * incX], test.want[i]))
            }
            checkValidIncGuard(xg, xGdVal, incX, gdln)
        }
    }
}

@Test
public func testScalIncTo() {
    const xGdVal = -1.0
    const dstGdVal = 0.5
    let gdln = 4
    let rng = Random(42)
    for (i in 0..scalTests.size) {
        let test = scalTests[i]
        let n = test.x.size
        for (inc in newIncSet(1, 2, 3, 4, 7, 10)) {
            let xg = guardIncVector(test.x, xGdVal, inc.x, gdln)
            let dg = guardIncVector(randSlice(n, 1, rng), dstGdVal, inc.y, gdln)
            let x = xg[gdln..xg.size - gdln]
            let dst = dg[gdln..dg.size - gdln]

            scalIncTo(dst, inc.y, test.alpha, x, n, inc.x)

            for (i in 0..test.want.size) {
                @Expect(same(dst[i * inc.y], test.want[i]))
            }
            checkValidIncGuard(xg, xGdVal, inc.x, gdln)
            checkValidIncGuard(dg, dstGdVal, inc.y, gdln)
            @Expect(equalStrided(test.x, x, inc.x))
        }
    }
}
