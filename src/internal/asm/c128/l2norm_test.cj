package cjnum.internal.asm.c128

import cjnum.complex.*
import std.math.*
import std.unittest.*
import std.unittest.testmacro.*

func naNWith(payload: UInt64): Complex128 {
    const nanBits: UInt64 = 0x7ff8000000000000
    const nanMask: UInt64 = 0xfff8000000000000
    var nan = Float64.fromBits(nanBits | (payload & (!nanMask)))
    return Complex128(nan, nan)
}

class L2NormUnitaryTestStruct {
    L2NormUnitaryTestStruct(
        let want!: Float64,
        var x!: Array<Complex128>
    ) {}
}

@Test
public func testL2NormUnitary() {
    const tol = 1e-15

    let src_gd: Complex128 = Complex128(1.0)
    let tests = [
        L2NormUnitaryTestStruct(want: 0.0, x: []),
        L2NormUnitaryTestStruct(want: 2.0, x: [Complex128(2.0)]),
        L2NormUnitaryTestStruct(want: 2.0, x: [Complex128(0.0, 2.0)]),
        L2NormUnitaryTestStruct(want: sqrt(8.0), x: [Complex128(2.0, 2.0)]),
        L2NormUnitaryTestStruct(want: 3.7416573867739413, x: [Complex128(1.0), Complex128(2.0), Complex128(3.0)]),
        L2NormUnitaryTestStruct(want: 3.7416573867739413, x: [Complex128(-1.0), Complex128(-2.0), Complex128(-3.0)]),
        L2NormUnitaryTestStruct(want: 3.7416573867739413,
            x: [Complex128(0.0, 1.0), Complex128(0.0, 2.0), Complex128(0.0, 3.0)]),
        L2NormUnitaryTestStruct(want: 3.7416573867739413,
            x: [Complex128(0.0, -1.0), Complex128(0.0, -2.0), Complex128(0.0, -3.0)]),
        L2NormUnitaryTestStruct(want: sqrt(28.0), x: [Complex128(1.0, 1.0), Complex128(2.0, 2.0), Complex128(3.0, 3.0)]),
        L2NormUnitaryTestStruct(
            want: sqrt(28.0),
            x: [Complex128(-1.0, -1.0), Complex128(-2.0, -2.0), Complex128(-3.0, -3.0)]
        ),
        L2NormUnitaryTestStruct(want: Float64.NaN, x: [Complex128.NaN]),
        L2NormUnitaryTestStruct(want: Float64.NaN,
            x: [Complex128(1.0, 0.0), Complex128(inf), Complex128(3.0, 0.0), naNWith(25), Complex128(5.0, 0.0)]),
        L2NormUnitaryTestStruct(want: 17.88854381999832,
            x: [Complex128(8.0), Complex128(-8.0), Complex128(8.0), Complex128(-8.0), Complex128(8.0)]),
        L2NormUnitaryTestStruct(
            want: 2.23606797749979,
            x: [Complex128(0.0), Complex128(1.0), Complex128(0.0), Complex128(-1.0), Complex128(0.0), Complex128(1.0),
                Complex128(0.0), Complex128(-1.0), Complex128(0.0), Complex128(1.0)]
        ),
        L2NormUnitaryTestStruct(
            want: 17.88854381999832,
            x: [Complex128(0.0, 8.0), Complex128(0.0, -8.0), Complex128(0.0, 8.0), Complex128(0.0, -8.0),
                Complex128(0.0, 8.0)]
        ),
        L2NormUnitaryTestStruct(
            want: 2.23606797749979,
            x: [Complex128(0.0), Complex128(0.0, 1.0), Complex128(0.0), Complex128(0.0, -1.0), Complex128(0.0),
                Complex128(0.0, 1.0), Complex128(0.0), Complex128(0.0, -1.0), Complex128(0.0), Complex128(0.0, 1.0)]
        ),
        L2NormUnitaryTestStruct(
            want: sqrt(640.0),
            x: [Complex128(8.0, 8.0), Complex128(-8.0, -8.0), Complex128(8.0, 8.0), Complex128(-8.0, -8.0),
                Complex128(8.0, 8.0)]
        ),
        L2NormUnitaryTestStruct(
            want: sqrt(10.0),
            x: [Complex128(0.0), Complex128(1.0, 1.0), Complex128(0.0), Complex128(-1.0, -1.0), Complex128(0.0),
                Complex128(1.0, 1.0), Complex128(0.0), Complex128(-1.0, -1.0), Complex128(0.0), Complex128(1.0, 1.0)]
        )
    ]
    for (j in 0..tests.size) {
        let g_ln = 4 + j % 2
        let x_ = guardVector(tests[j].x, src_gd, g_ln)
        let src = x_[g_ln..x_.size - g_ln]
        let ret = l2NormUnitary(src)
        @Expect(sameApprox(ret, tests[j].want, tol))
        @Expect(isValidGuard(x_, src_gd, g_ln))
    }
}

class L2DistanceUnitaryTestStruct {
    L2DistanceUnitaryTestStruct(
        let want!: Float64,
        var x!: Array<Complex128>,
        var y!: Array<Complex128>
    ) {}
}

@Test
public func testL2DistanceUnitary() {
    const tol = 1e-15

    let src_gd: Complex128 = Complex128(1.0)
    let tests = [
        L2DistanceUnitaryTestStruct(want: 0.0, x: [], y: []),
        L2DistanceUnitaryTestStruct(want: 2.0, x: [Complex128(3.0, 0.0)], y: [Complex128(1.0, 0.0)]),
        L2DistanceUnitaryTestStruct(want: 2.0, x: [Complex128(0.0, 3.0)], y: [Complex128(0.0, 1.0)]),
        L2DistanceUnitaryTestStruct(
            want: 3.7416573867739413,
            x: [Complex128(2.0, 0.0), Complex128(4.0, 0.0), Complex128(6.0, 0.0)],
            y: [Complex128(1.0, 0.0), Complex128(2.0, 0.0), Complex128(3.0, 0.0)]
        ),
        L2DistanceUnitaryTestStruct(
            want: 3.7416573867739413,
            x: [Complex128(1.0, 0.0), Complex128(2.0, 0.0), Complex128(3.0, 0.0)],
            y: [Complex128(2.0, 0.0), Complex128(4.0, 0.0), Complex128(6.0, 0.0)]
        ),
        L2DistanceUnitaryTestStruct(
            want: 3.7416573867739413,
            x: [Complex128(0.0, 2.0), Complex128(0.0, 4.0), Complex128(0.0, 6.0)],
            y: [Complex128(0.0, 1.0), Complex128(0.0, 2.0), Complex128(0.0, 3.0)]
        ),
        L2DistanceUnitaryTestStruct(
            want: 3.7416573867739413,
            x: [Complex128(0.0, 1.0), Complex128(0.0, 2.0), Complex128(0.0, 3.0)],
            y: [Complex128(0.0, 2.0), Complex128(0.0, 4.0), Complex128(0.0, 6.0)]
        ),
        L2DistanceUnitaryTestStruct(
            want: sqrt(28.0),
            x: [Complex128(2.0, 2.0), Complex128(4.0, 4.0), Complex128(6.0, 6.0)],
            y: [Complex128(1.0, 1.0), Complex128(2.0, 2.0), Complex128(3.0, 3.0)]
        ),
        L2DistanceUnitaryTestStruct(
            want: sqrt(28.0),
            x: [Complex128(1.0, 1.0), Complex128(2.0, 2.0), Complex128(3.0, 3.0)],
            y: [Complex128(2.0, 2.0), Complex128(4.0, 4.0), Complex128(6.0, 6.0)]
        ),
        L2DistanceUnitaryTestStruct(want: Float64.NaN, x: [Complex128.NaN], y: [Complex128(0.0, 0.0)]),
        L2DistanceUnitaryTestStruct(
            want: 17.88854381999832,
            x: [Complex128(9.0, 0.0), Complex128(-9.0, 0.0), Complex128(9.0, 0.0), Complex128(-9.0, 0.0),
                Complex128(9.0, 0.0)],
            y: [Complex128(1.0, 0.0), Complex128(-1.0, 0.0), Complex128(1.0, 0.0), Complex128(-1.0, 0.0),
                Complex128(1.0, 0.0)]
        ),
        L2DistanceUnitaryTestStruct(
            want: 2.23606797749979,
            x: [Complex128(0.0, 0.0), Complex128(1.0, 0.0), Complex128(0.0, 0.0), Complex128(-1.0, 0.0),
                Complex128(0.0, 0.0), Complex128(1.0, 0.0), Complex128(0.0, 0.0), Complex128(-1.0, 0.0),
                Complex128(0.0, 0.0), Complex128(1.0, 0.0)],
            y: [Complex128(0.0, 0.0), Complex128(2.0, 0.0), Complex128(0.0, 0.0), Complex128(-2.0, 0.0),
                Complex128(0.0, 0.0), Complex128(2.0, 0.0), Complex128(0.0, 0.0), Complex128(-2.0, 0.0),
                Complex128(0.0, 0.0), Complex128(2.0, 0.0)]
        ),
        L2DistanceUnitaryTestStruct(
            want: 17.88854381999832,
            x: [Complex128(0.0, 9.0), Complex128(0.0, -9.0), Complex128(0.0, 9.0), Complex128(0.0, -9.0),
                Complex128(0.0, 9.0)],
            y: [Complex128(0.0, 1.0), Complex128(0.0, -1.0), Complex128(0.0, 1.0), Complex128(0.0, -1.0),
                Complex128(0.0, 1.0)]
        ),
        L2DistanceUnitaryTestStruct(
            want: 2.23606797749979,
            x: [Complex128(0.0, 0.0), Complex128(0.0, 1.0), Complex128(0.0, 0.0), Complex128(0.0, -1.0),
                Complex128(0.0, 0.0), Complex128(0.0, 1.0), Complex128(0.0, 0.0), Complex128(0.0, -1.0),
                Complex128(0.0, 0.0), Complex128(0.0, 1.0)],
            y: [Complex128(0.0, 0.0), Complex128(0.0, 2.0), Complex128(0.0, 0.0), Complex128(0.0, -2.0),
                Complex128(0.0, 0.0), Complex128(0.0, 2.0), Complex128(0.0, 0.0), Complex128(0.0, -2.0),
                Complex128(0.0, 0.0), Complex128(0.0, 2.0)]
        ),
        L2DistanceUnitaryTestStruct(
            want: sqrt(640.0),
            x: [Complex128(9.0, 9.0), Complex128(-9.0, -9.0), Complex128(9.0, 9.0), Complex128(-9.0, -9.0),
                Complex128(9.0, 9.0)],
            y: [Complex128(1.0, 1.0), Complex128(-1.0, -1.0), Complex128(1.0, 1.0), Complex128(-1.0, -1.0),
                Complex128(1.0, 1.0)]
        ),
        L2DistanceUnitaryTestStruct(
            want: sqrt(10.0),
            x: [Complex128(0.0, 0.0), Complex128(1.0, 1.0), Complex128(0.0, 0.0), Complex128(-1.0, -1.0),
                Complex128(0.0, 0.0), Complex128(1.0, 1.0), Complex128(0.0, 0.0), Complex128(-1.0, -1.0),
                Complex128(0.0, 0.0), Complex128(1.0, 1.0)],
            y: [Complex128(0.0, 0.0), Complex128(2.0, 2.0), Complex128(0.0, 0.0), Complex128(-2.0, -2.0),
                Complex128(0.0, 0.0), Complex128(2.0, 2.0), Complex128(0.0, 0.0), Complex128(-2.0, -2.0),
                Complex128(0.0, 0.0), Complex128(2.0, 2.0)]
        )
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

func netlib(x: Array<Complex128>): Float64 {
    var scale: Float64 = 0.0
    var sumSquares: Float64 = 1.0
    for (v in x) {
        if (v == Complex128(0.0)) {
            continue
        }
        let absxi = c128Abs(v)
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
    func xSlice(ln: Int64): Array<Complex128> {
        let x: Array<Complex128> = Array(ln) {i => Complex128(Float64(i))}
        x[0] = Complex128(4.0)
        return x[..ln]
    }

    @Bench[x_ in xSlice]
    @Configure[]
    func testL2NormUnitaryNetlib(x_: Array<Complex128>): Unit {
        netlib(x_)
    }

    @Bench[x_ in xSlice]
    @Configure[]
    func testL2NormUnitary(x_: Array<Complex128>): Unit {
        l2NormUnitary(x_)
    }
}
