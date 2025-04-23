package cjnum.internal.asm.f64

import std.math.*
import std.unittest.*
import std.unittest.testmacro.*

const nanBits: UInt64 = 0x7ff8000000000000
const nanMask: UInt64 = 0xfff8000000000000

public func naNWith(payload: UInt64): Float64 {
    return Float64.fromBits(nanBits | (payload & (!nanMask)))
}

class L2NormUnitaryTestStruct {
    L2NormUnitaryTestStruct(
        let want!: Float64,
        var x!: Array<Float64>
    ) {}
}

@Test
public func testL2NormUnitary() {
    const tol = 1e-15

    let src_gd: Float64 = 1.0
    let tests = [
        L2NormUnitaryTestStruct(want: 0.0, x: []),
        L2NormUnitaryTestStruct(want: 2.0, x: [2.0]),
        L2NormUnitaryTestStruct(want: 3.7416573867739413, x: [1.0, 2.0, 3.0]),
        L2NormUnitaryTestStruct(want: 3.7416573867739413, x: [-1.0, -2.0, -3.0]),
        L2NormUnitaryTestStruct(want: nan, x: [nan]),
        L2NormUnitaryTestStruct(want: nan, x: [1.0, inf, 3.0, naNWith(25), 5.0]),
        L2NormUnitaryTestStruct(want: 17.88854381999832, x: [8.0, -8.0, 8.0, -8.0, 8.0]),
        L2NormUnitaryTestStruct(want: 2.23606797749979, x: [0.0, 1.0, 0.0, -1.0, 0.0, 1.0, 0.0, -1.0, 0.0, 1.0])
    ]
    for (j in 0..tests.size) {
        let g_ln = 4 + j % 2
        tests[j].x = guardVector(tests[j].x, src_gd, g_ln)
        let src = tests[j].x[g_ln..tests[j].x.size - g_ln]
        let ret = l2NormUnitary(src)
        @Expect(sameApprox(ret, tests[j].want, tol))
        @Expect(isValidGuard(tests[j].x, src_gd, g_ln))
    }
}

class L2NormIncTestStruct {
    L2NormIncTestStruct(
        let inc!: Int64,
        let want!: Float64,
        var x!: Array<Float64>
    ) {}
}

@Test
public func testL2NormInc() {
    const tol = 1e-15

    let src_gd: Float64 = 1.0
    let tests = [
        L2NormIncTestStruct(inc: 2, want: 0.0, x: []),
        L2NormIncTestStruct(inc: 3, want: 2.0, x: [2.0]),
        L2NormIncTestStruct(inc: 10, want: 3.7416573867739413, x: [1.0, 2.0, 3.0]),
        L2NormIncTestStruct(inc: 5, want: 3.7416573867739413, x: [-1.0, -2.0, -3.0]),
        L2NormIncTestStruct(inc: 3, want: nan, x: [nan]),
        L2NormIncTestStruct(inc: 15, want: 17.88854381999832, x: [8.0, -8.0, 8.0, -8.0, 8.0]),
        L2NormIncTestStruct(inc: 1, want: 2.23606797749979, x: [0.0, 1.0, 0.0, -1.0, 0.0, 1.0, 0.0, -1.0, 0.0, 1.0])
    ]
    for (j in 0..tests.size) {
        let g_ln = 4 + j % 2
        let ln = tests[j].x.size
        tests[j].x = guardIncVector(tests[j].x, src_gd, tests[j].inc, g_ln)
        let src = tests[j].x[g_ln..tests[j].x.size - g_ln]
        let ret = l2NormInc(src, ln, tests[j].inc)
        @Expect(sameApprox(ret, tests[j].want, tol))
        @Expect(checkValidIncGuard(tests[j].x, src_gd, tests[j].inc, g_ln))
    }
}

class L2DistanceUnitaryTestStruct {
    L2DistanceUnitaryTestStruct(
        let want!: Float64,
        var x!: Array<Float64>,
        var y!: Array<Float64>
    ) {}
}

@Test
public func testL2DistanceUnitary() {
    const tol = 1e-15

    let src_gd: Float64 = 1.0
    let tests = [
        L2DistanceUnitaryTestStruct(want: 0.0, x: [], y: []),
        L2DistanceUnitaryTestStruct(want: 2.0, x: [3.0], y: [1.0]),
        L2DistanceUnitaryTestStruct(want: 3.7416573867739413, x: [2.0, 4.0, 6.0], y: [1.0, 2.0, 3.0]),
        L2DistanceUnitaryTestStruct(want: 3.7416573867739413, x: [1.0, 2.0, 3.0], y: [2.0, 4.0, 6.0]),
        L2DistanceUnitaryTestStruct(want: nan, x: [nan], y: [0.0]),
        L2DistanceUnitaryTestStruct(want: 17.88854381999832, x: [9.0, -9.0, 9.0, -9.0, 9.0], y: [1.0, -1.0, 1.0, -1.0,
            1.0]),
        L2DistanceUnitaryTestStruct(want: 2.23606797749979, x: [0.0, 1.0, 0.0, -1.0, 0.0, 1.0, 0.0, -1.0, 0.0, 1.0],
            y: [0.0, 2.0, 0.0, -2.0, 0.0, 2.0, 0.0, -2.0, 0.0, 2.0])
    ]
    for (j in 0..tests.size) {
        let g_ln = 4 + j % 2
        tests[j].x = guardVector(tests[j].x, src_gd, g_ln)
        tests[j].y = guardVector(tests[j].y, src_gd, g_ln)
        let srcX = tests[j].x[g_ln..tests[j].x.size - g_ln]
        let srcY = tests[j].y[g_ln..tests[j].y.size - g_ln]
        let ret = l2DistanceUnitary(srcX, srcY)
        @Expect(sameApprox(ret, tests[j].want, tol))
        @Expect(isValidGuard(tests[j].x, src_gd, g_ln))
    }
}

func netlib(x: Array<Float64>): Float64 {
    var scale: Float64 = 0.0
    var sumSquares: Float64 = 1.0
    for (v in x) {
        if (v == 0.0) {
            continue
        }
        let absxi = abs(v)
        if (absxi.isNaN()) {
            return Float64.NaN
        }
        if (scale < absxi) {
            let s = scale / absxi
            sumSquares = 1.0 + sumSquares * s * s
            scale = absxi
        } else {
            let s = absxi / scale
            sumSquares += s * s
        }
    }
    if (scale == Float64.Inf) {
        return Float64.Inf
    }
    return scale * sqrt(sumSquares)
}

@Test
public class BenchmarkL2NormNetlibTestClass {
    let benchmarkSizes = [1, 3, 10, 30, Int64(1e2), Int64(3e2), Int64(1e3), Int64(3e3), Int64(1e4), Int64(3e4),
        Int64(1e5)]
    
    @Strategy[ln in benchmarkSizes]
    func xSlice(ln: Int64): Array<Float64> {
        let num = randomSlice(1, 1)[0]
        var x = Array(ln) { i => Float64(i) }
        x[0] = num
        return x[..ln]
    }

    @Bench[x_ in xSlice]
    @Configure[]
    func testL2NormUnitaryNetlib(x_: Array<Float64>): Unit {
        netlib(x_)
    }

    @Bench[x_ in xSlice]
    @Configure[]
    func testL2NormUnitary(x_: Array<Float64>): Unit {
        l2NormUnitary(x_)
    }
}
