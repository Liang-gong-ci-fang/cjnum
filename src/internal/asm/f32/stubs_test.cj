package cjnum.internal.asm.f32

import std.unittest.*
import std.unittest.testmacro.*

class testStruct {
    testStruct(
        let incX: Int64,
        let incY: Int64,
        let incDst: Int64,
        let ix: Int64,
        let iy: Int64,
        let idst: Int64,
        let a: Float32,
        var dst: Array<Float32>,
        var x: Array<Float32>,
        var y: Array<Float32>,
        let ex: Array<Float32>
    ) {}
}

@Test
public func testAxpyUnitary() {
    let tests = [
        testStruct(2, 2, 3, 0, 0, 0, 3.0, [5.0], [2.0], [1.0], [7.0]),
        testStruct(2, 2, 3, 0, 0, 0, 5.0, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [1.0, 1.0, 1.0], [1.0, 1.0, 1.0]),
        testStruct(2, 2, 3, 0, 0, 0, 5.0, [0.0, 0.0, 0.0], [0.0, 0.0], [1.0, 1.0, 1.0], [1.0, 1.0]),
        testStruct(2, 2, 3, 0, 0, 0, -1.0, [-1.0, -1.0, -1.0], [1.0, 1.0, 1.0], [1.0, 2.0, 1.0], [0.0, 1.0, 0.0]),
        testStruct(2, 2, 3, 0, 0, 0, -1.0, [1.0, 1.0, 1.0], [1.0, 2.0, 1.0], [-1.0, -2.0, -1.0], [-2.0, -4.0, -2.0]),
        testStruct(2, 2, 3, 0, 0, 0, 2.5, [1.0, 1.0, 1.0, 1.0, 1.0], [1.0, 2.0, 3.0, 2.0, 1.0], [0.0, 0.0, 0.0, 0.0, 0.0],
            [2.5, 5.0, 7.5, 5.0, 2.5]),
        testStruct(2, 2, 3, 0, 0, 0, 16.5, Array<Float32>(20, item: 0.0),
            [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5],
            [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0],
            [9.25, 10.25, 11.25, 12.25, 13.25, 14.25, 15.25, 16.25, 17.25, 18.25, 9.25, 10.25, 11.25, 12.25, 13.25,
            14.25, 15.25, 16.25, 17.25, 18.25]),
        testStruct(2, 2, 3, 0, 0, 0, 16.5, Array<Float32>(10, item: 0.0), [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
            0.5], [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0],
            [9.25, 10.25, 11.25, 12.25, 13.25, 14.25, 15.25, 16.25, 17.25, 18.25])
    ]
    for (j in 0..tests.size) {
        let v = tests[j]
        let gdLn = 4 + j % 2
        v.x = guardVector(v.x, 1.0, gdLn)
        v.y = guardVector(v.y, 1.0, gdLn)
        let x = v.x[gdLn..v.x.size - gdLn]
        let y = v.y[gdLn..v.y.size - gdLn]
        axpyUnitary(v.a, x, y)
        for (i in 0..x.size) {
            @Expect(same(y[i], v.ex[i]))
        }
        @Expect(isValidGuard(v.x, 1.0, gdLn))
        @Expect(isValidGuard(v.y, 1.0, gdLn))
    }
}

@Test
public func testAxpyUnitaryTo() {
    let tests = [
        testStruct(2, 2, 3, 0, 0, 0, 3.0, [5.0], [2.0], [1.0], [7.0]),
        testStruct(2, 2, 3, 0, 0, 0, 5.0, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [1.0, 1.0, 1.0], [1.0, 1.0, 1.0]),
        testStruct(2, 2, 3, 0, 0, 0, 5.0, [0.0, 0.0, 0.0], [0.0, 0.0], [1.0, 1.0, 1.0], [1.0, 1.0]),
        testStruct(2, 2, 3, 0, 0, 0, -1.0, [-1.0, -1.0, -1.0], [1.0, 1.0, 1.0], [1.0, 2.0, 1.0], [0.0, 1.0, 0.0]),
        testStruct(2, 2, 3, 0, 0, 0, -1.0, [1.0, 1.0, 1.0], [1.0, 2.0, 1.0], [-1.0, -2.0, -1.0], [-2.0, -4.0, -2.0]),
        testStruct(2, 2, 3, 0, 0, 0, 2.5, [1.0, 1.0, 1.0, 1.0, 1.0], [1.0, 2.0, 3.0, 2.0, 1.0], [0.0, 0.0, 0.0, 0.0, 0.0],
            [2.5, 5.0, 7.5, 5.0, 2.5]),
        testStruct(2, 2, 3, 0, 0, 0, 16.5, Array<Float32>(20, item: 0.0),
            [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5],
            [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0],
            [9.25, 10.25, 11.25, 12.25, 13.25, 14.25, 15.25, 16.25, 17.25, 18.25, 9.25, 10.25, 11.25, 12.25, 13.25,
            14.25, 15.25, 16.25, 17.25, 18.25]),
        testStruct(2, 2, 3, 0, 0, 0, 16.5, Array<Float32>(10, item: 0.0), [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
            0.5], [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0],
            [9.25, 10.25, 11.25, 12.25, 13.25, 14.25, 15.25, 16.25, 17.25, 18.25])
    ]
    for (j in 0..tests.size) {
        let v = tests[j]
        let gdLn = 4 + j % 2
        v.x = guardVector(v.x, 1.0, gdLn)
        v.y = guardVector(v.y, 1.0, gdLn)
        v.dst = guardVector(v.dst, 0.0, gdLn)
        let x = v.x[gdLn..v.x.size - gdLn]
        let y = v.y[gdLn..v.y.size - gdLn]
        let dst = v.dst[gdLn..v.dst.size - gdLn]
        axpyUnitaryTo(dst, v.a, x, y)
        for (i in 0..x.size) {
            @Expect(same(dst[i], v.ex[i]))
        }
        @Expect(isValidGuard(v.x, 1.0, gdLn))
        @Expect(isValidGuard(v.y, 1.0, gdLn))
        @Expect(isValidGuard(v.dst, 0.0, gdLn))
    }
}

@Test
public func testAxpyInc() {
    let tests = [
        testStruct(2, 2, 3, 0, 0, 0, 3.0, [5.0], [2.0], [1.0], [7.0]),
        testStruct(2, 2, 3, 0, 0, 0, 5.0, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [1.0, 1.0, 1.0], [1.0, 1.0, 1.0]),
        testStruct(2, 2, 3, 0, 0, 0, 5.0, [0.0, 0.0, 0.0], [0.0, 0.0], [1.0, 1.0, 1.0], [1.0, 1.0]),
        testStruct(2, 2, 3, 0, 0, 0, -1.0, [-1.0, -1.0, -1.0], [1.0, 1.0, 1.0], [1.0, 2.0, 1.0], [0.0, 1.0, 0.0]),
        testStruct(2, 2, 3, 0, 0, 0, -1.0, [1.0, 1.0, 1.0], [1.0, 2.0, 1.0], [-1.0, -2.0, -1.0], [-2.0, -4.0, -2.0]),
        testStruct(2, 2, 3, 0, 0, 0, 2.5, [1.0, 1.0, 1.0, 1.0, 1.0], [1.0, 2.0, 3.0, 2.0, 1.0], [0.0, 0.0, 0.0, 0.0, 0.0],
            [2.5, 5.0, 7.5, 5.0, 2.5]),
        testStruct(2, 2, 3, 0, 0, 0, 16.5, Array<Float32>(20, item: 0.0),
            [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5],
            [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0],
            [9.25, 10.25, 11.25, 12.25, 13.25, 14.25, 15.25, 16.25, 17.25, 18.25, 9.25, 10.25, 11.25, 12.25, 13.25,
            14.25, 15.25, 16.25, 17.25, 18.25]),
        testStruct(2, 2, 3, 0, 0, 0, 16.5, Array<Float32>(10, item: 0.0), [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
            0.5], [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0],
            [9.25, 10.25, 11.25, 12.25, 13.25, 14.25, 15.25, 16.25, 17.25, 18.25])
    ]
    for (j in 0..tests.size) {
        let v = tests[j]
        let gdLn = 4 + j % 2
        v.x = guardIncVector(v.x, 1.0, v.incX, gdLn)
        v.y = guardIncVector(v.y, 1.0, v.incY, gdLn)
        let x = v.x[gdLn..v.x.size - gdLn]
        let y = v.y[gdLn..v.y.size - gdLn]
        axpyInc(v.a, x, y, v.ex.size, v.incX, v.incY, v.ix, v.iy)
        for (i in 0..v.ex.size) {
            @Expect(same(y[i*v.incY], v.ex[i]))
        }
        @Expect(checkValidIncGuard(v.x, 1.0, v.incX, gdLn))
        @Expect(checkValidIncGuard(v.y, 1.0, v.incY, gdLn))
    }
}

@Test
public func testAxpyIncTo() {
    let tests = [
        testStruct(2, 2, 3, 0, 0, 0, 3.0, [5.0], [2.0], [1.0], [7.0]),
        testStruct(2, 2, 3, 0, 0, 0, 5.0, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [1.0, 1.0, 1.0], [1.0, 1.0, 1.0]),
        testStruct(2, 2, 3, 0, 0, 0, 5.0, [0.0, 0.0, 0.0], [0.0, 0.0], [1.0, 1.0, 1.0], [1.0, 1.0]),
        testStruct(2, 2, 3, 0, 0, 0, -1.0, [-1.0, -1.0, -1.0], [1.0, 1.0, 1.0], [1.0, 2.0, 1.0], [0.0, 1.0, 0.0]),
        testStruct(2, 2, 3, 0, 0, 0, -1.0, [1.0, 1.0, 1.0], [1.0, 2.0, 1.0], [-1.0, -2.0, -1.0], [-2.0, -4.0, -2.0]),
        testStruct(2, 2, 3, 0, 0, 0, 2.5, [1.0, 1.0, 1.0, 1.0, 1.0], [1.0, 2.0, 3.0, 2.0, 1.0], [0.0, 0.0, 0.0, 0.0, 0.0],
            [2.5, 5.0, 7.5, 5.0, 2.5]),
        testStruct(2, 2, 3, 0, 0, 0, 16.5, Array<Float32>(20, item: 0.0),
            [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5],
            [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0],
            [9.25, 10.25, 11.25, 12.25, 13.25, 14.25, 15.25, 16.25, 17.25, 18.25, 9.25, 10.25, 11.25, 12.25, 13.25,
            14.25, 15.25, 16.25, 17.25, 18.25]),
        testStruct(2, 2, 3, 0, 0, 0, 16.5, Array<Float32>(10, item: 0.0), [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
            0.5], [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0],
            [9.25, 10.25, 11.25, 12.25, 13.25, 14.25, 15.25, 16.25, 17.25, 18.25])
    ]
    for (j in 0..tests.size) {
        let v = tests[j]
        let gdLn = 4 + j % 2
        v.x = guardIncVector(v.x, 1.0, v.incX, gdLn)
        v.y = guardIncVector(v.y, 1.0, v.incY, gdLn)
        v.dst = guardIncVector(v.dst, 0.0, v.incDst, gdLn)
        let x = v.x[gdLn..v.x.size - gdLn]
        let y = v.y[gdLn..v.y.size - gdLn]
        let dst = v.dst[gdLn..v.dst.size - gdLn]
        axpyIncTo(dst, v.incDst, v.idst, v.a, x, y, v.ex.size, v.incX, v.incY, v.ix, v.iy)
        for (i in 0..v.ex.size) {
            @Expect(same(dst[i*v.incDst], v.ex[i]))
        }
        @Expect(checkValidIncGuard(v.x, 1.0, v.incX, gdLn))
        @Expect(checkValidIncGuard(v.y, 1.0, v.incY, gdLn))
        @Expect(checkValidIncGuard(v.dst, 0.0, v.incDst, gdLn))
    }
}

class sumTestStruct {
    sumTestStruct(
        let src: Array<Float32>,
        let expect: Float32
    ) {}
}

@Test
public func testSum() {
    let srcGd: Float32 = -1.0
    let sumTests = [
        sumTestStruct([], 0.0),
        sumTestStruct([1.0], 1.0),
        sumTestStruct([nan], nan),
        sumTestStruct([1.0, 2.0, 3.0], 6.0),
        sumTestStruct([1.0, -4.0, 3.0], 0.0),
        sumTestStruct([1.0, 2.0, 3.0, 4.0], 10.0),
        sumTestStruct([1.0, 1.0, nan, 1.0, 1.0], nan),
        sumTestStruct([inf, 4.0, nan, -inf, 9.0], nan),
        sumTestStruct([1.0, 1.0, 1.0, 1.0, 9.0, 1.0, 1.0, 1.0, 2.0, 1.0, 1.0, 1.0, 1.0, 1.0, 5.0, 1.0], 29.0),
        sumTestStruct(
            [1.0, 1.0, 1.0, 1.0, 9.0, 1.0, 1.0, 1.0, 2.0, 1.0, 1.0, 1.0, 1.0, 1.0, 5.0, 11.0, 1.0, 1.0, 1.0, 9.0, 1.0,
                1.0, 1.0, 2.0, 1.0, 1.0, 1.0, 1.0, 1.0, 5.0, 1.0],
            67.0
        )
    ]
    for (j in 0..sumTests.size) {
        let v = sumTests[j]
        for (i in 0..4) {
            let gdLn = 4 + j % 4 + i
            let gsrc = guardVector(v.src, srcGd, gdLn)
            let src = gsrc[gdLn..gsrc.size - gdLn]
            let ret = sum(src)
            @Expect(same(ret, v.expect))
            @Expect(isValidGuard(gsrc, srcGd, gdLn))
        }
    }
}
