package cjnum.internal.asm.f64

import cjnum.floats.scalar.*
import std.math.*
import std.random.*
import std.unittest.*
import std.unittest.testmacro.*

class DotUnitaryTestStruct {
    DotUnitaryTestStruct(
        let xData!: Array<Float64>,
        let yData!: Array<Float64>,
        let want!: Float64
    ) {}
}

@Test
public func testDotUnitary() {
    let tests = [
        DotUnitaryTestStruct(xData: [2.0], yData: [-3.0], want: -6.0),
        DotUnitaryTestStruct(xData: [2.0, 3.0], yData: [-3.0, 4.0], want: 6.0),
        DotUnitaryTestStruct(xData: [2.0, 3.0, -4.0], yData: [-3.0, 4.0, 5.0], want: -14.0),
        DotUnitaryTestStruct(xData: [2.0, 3.0, -4.0, -5.0], yData: [-3.0, 4.0, 5.0, -6.0], want: 16.0),
        DotUnitaryTestStruct(xData: [0.0, 2.0, 3.0, -4.0, -5.0], yData: [0.0, -3.0, 4.0, 5.0, -6.0], want: 16.0),
        DotUnitaryTestStruct(xData: [0.0, 0.0, 2.0, 3.0, -4.0, -5.0], yData: [0.0, 1.0, -3.0, 4.0, 5.0, -6.0],
            want: 16.0),
        DotUnitaryTestStruct(xData: [0.0, 0.0, 1.0, 1.0, 2.0, -3.0, -4.0], yData: [0.0, 1.0, 0.0, 3.0, -4.0, 5.0, -6.0],
            want: 4.0),
        DotUnitaryTestStruct(xData: [0.0, 0.0, 1.0, 1.0, 2.0, -3.0, -4.0, 5.0],
            yData: [0.0, 1.0, 0.0, 3.0, -4.0, 5.0, -6.0, 7.0], want: 39.0)
    ]
    for (test in tests) {
        let xRestlt = newGuardedVector(test.xData, 1)
        let x = xRestlt[0]
        let xFront = xRestlt[1]
        let xBack = xRestlt[2]
        let yRestlt = newGuardedVector(test.yData, 1)
        let y = yRestlt[0]
        let yFront = yRestlt[1]
        let yBack = yRestlt[2]
        let got = dotUnitary(x, y)
        @Expect(allNaN(xFront) || allNaN(xBack))
        @Expect(allNaN(yFront) || allNaN(yBack))
        @Expect(equalStrided(test.xData, x, 1))
        @Expect(equalStrided(test.yData, y, 1))
        @Expect(!same(got,Float64.NaN))
        @Expect(got,test.want)
    }
}

class DotIncTestStruct {
    DotIncTestStruct(
        let xData!: Array<Float64>,
        let yData!: Array<Float64>,
        let want!: Float64,
        let wantRev!: Float64
    ) {}
}

@Test
public func testDotInc() {
    let tests = [
        DotIncTestStruct(xData: [2.0], yData: [-3.0], want: -6.0, wantRev: -6.0),
        DotIncTestStruct(xData: [2.0, 3.0], yData: [-3.0, 4.0], want: 6.0, wantRev: -1.0),
        DotIncTestStruct(xData: [2.0, 3.0, -4.0], yData: [-3.0, 4.0, 5.0], want: -14.0, wantRev: 34.0),
        DotIncTestStruct(xData: [2.0, 3.0, -4.0, -5.0], yData: [-3.0, 4.0, 5.0, -6.0], want: 16.0, wantRev: 2.0),
        DotIncTestStruct(xData: [0.0, 2.0, 3.0, -4.0, -5.0], yData: [0.0, -3.0, 4.0, 5.0, -6.0], want: 16.0,
            wantRev: 34.0),
        DotIncTestStruct(xData: [0.0, 0.0, 2.0, 3.0, -4.0, -5.0], yData: [0.0, 1.0, -3.0, 4.0, 5.0, -6.0], want: 16.0,
            wantRev: -5.0),
        DotIncTestStruct(xData: [0.0, 0.0, 1.0, 1.0, 2.0, -3.0, -4.0], yData: [0.0, 1.0, 0.0, 3.0, -4.0, 5.0, -6.0],
            want: 4.0, wantRev: -4.0),
        DotIncTestStruct(xData: [0.0, 0.0, 1.0, 1.0, 2.0, -3.0, -4.0, 5.0], yData: [0.0, 1.0, 0.0, 3.0, -4.0, 5.0, -6.0,
            7.0], want: 39.0, wantRev: 3.0)
    ]
    let range = [-7, -3, -2, -1, 1, 2, 3, 7]
    for (test in tests) {
        for (incX in range) {
            for (incY in range) {
                let n = test.xData.size
                let xRestlt = newGuardedVector(test.xData, incX)
                let x = xRestlt[0]
                let xFront = xRestlt[1]
                let xBack = xRestlt[2]
                let yRestlt = newGuardedVector(test.yData, incY)
                let y = yRestlt[0]
                let yFront = yRestlt[1]
                let yBack = yRestlt[2]
                var ix = 0
                var iy = 0
                if (incX < 0) {
                    ix = (-n + 1) * incX
                }
                if (incY < 0) {
                    iy = (-n + 1) * incY
                }
                let got = dotInc(x, y, n, incX, incY, ix, iy)
                
                @Expect(allNaN(xFront) || allNaN(xBack))
                @Expect(allNaN(yFront) || allNaN(yBack))
                @Expect(!nonStridedWrite(x, incX) || equalStrided(test.xData, x, incX))
                @Expect(!nonStridedWrite(y, incY) || equalStrided(test.yData, y, incY))
                @Expect(!same(got,Float64.NaN))
                var want = test.want
                if (incX * incY < 0) {
                    want = test.wantRev
                }
                @Expect(got,want)
            }
        }
    }
}

@Test
public class BenchmarkDotUnitaryN1TestClass {
    let benchmarkSizes = [1, 2, 3, 4, 10, Int64(1e2), Int64(1e3), Int64(1e4), Int64(1e5)]

    @Strategy[ln in benchmarkSizes]
    func dotSlice(ln: Int64): (Array<Float64>, Array<Float64>) {
        let x = Array<Float64>(ln, repeat: 0.0)
        let random = Random()
        for (i in 0..ln) {
            x[i] = random.nextFloat64()
        }
        let y = Array<Float64>(ln, repeat: 0.0)
        for (i in 0..ln) {
            y[i] = random.nextFloat64()
        }
        return (x[..ln], y[..ln])
    }

    @Bench[tests in dotSlice]
    @Configure[]
    func testBenchmarkDotUnitary(tests: (Array<Float64>, Array<Float64>)): Unit {
        dotUnitary(tests[0], tests[1])
    }

    let dotIncbenchmark = [[1, 1], [2, 1], [2, 2], [2, 4], [2, 10], [3, 1], [3, 2], [3, 4], [3, 10], [4, 1], [4, 2],
        [4, 4], [4, 10], [10, 1], [10, 2], [10, 4], [10, 10], [Int64(1e3), 1], [Int64(1e3), 2], [Int64(1e3), 4],
        [Int64(1e3), 10], [Int64(1e5), 1], [Int64(1e5), 2], [Int64(1e5), 4], [Int64(1e5), 10], [Int64(1e5), -1],
        [Int64(1e5), -2], [Int64(1e5), -4], [Int64(1e5), -10]]

    @Strategy[tests in dotIncbenchmark]
    func incSlice(tests: Array<Int64>): (Array<Float64>, Array<Float64>, Int64, Int64, Int64, Int64, Int64) {
        let n = tests[0]
        let inc = tests[1]
        let absInc = abs(inc)
        let ln = (n - 1) * absInc + 1
        let x = Array<Float64>(ln, repeat: 0.0)
        let y = Array<Float64>(ln, repeat: 0.0)
        var ini: Int64 = 0
        let random = Random()
        for (i in 0..x.size) {
            x[i] = random.nextFloat64()
        }
        for (i in 0..y.size) {
            y[i] = random.nextFloat64()
        }
        if (inc < 0) {
            ini = (-n + 1) * inc
        }
        return (x[..ln], y[..ln], n, inc, inc, ini, ini)
    }

    @Bench[incResult in incSlice]
    @Configure[]
    func testBenchmarkDotInc(incResult: (Array<Float64>, Array<Float64>, Int64, Int64, Int64, Int64, Int64)): Unit {
        dotInc(incResult[0], incResult[1], incResult[2], incResult[3], incResult[4], incResult[5], incResult[6])
    }
}
