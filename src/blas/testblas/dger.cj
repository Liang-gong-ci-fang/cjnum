package cjnum.blas.testblas

import std.math.*
import cjnum.blas.*
import cjnum.floats.*

class dgerTestStruct {
    dgerTestStruct(
        let name!: String,
        let a!: Array<Array<Float64>>,
        let m!: Int64,
        let n!: Int64,
        let x!: Array<Float64>,
        let y!: Array<Float64>,
        let incX!: Int64,
        let incY!: Int64,
        let want!: Array<Array<Float64>>
    ) {}
}

public func dgerTest(nFloat64: NFloat64): Bool {
    let tests = [
        dgerTestStruct(
            name: "M gt N inc 1",
            m: 5,
            n: 3,
            a: [
                [1.3, 2.4, 3.5],
                [2.6, 2.8, 3.3],
                [-1.3, -4.3, -9.7],
                [8.0, 9.0, -10.0],
                [-12.0, -14.0, -6.0]
            ],
            x: [-2.0, -3.0, 0.0, 1.0, 2.0],
            y: [-1.1, 5.0, 0.0],
            incX: 1,
            incY: 1,
            want: [
                [3.5, -7.6, 3.5], 
                [5.9, -12.2, 3.3], 
                [-1.3, -4.3, -9.7], 
                [6.9, 14.0, -10.0], 
                [-14.2, -4.0, -6.0]
            ]
        ),
        dgerTestStruct(
            name: "M eq N inc 1",
            m: 3,
            n: 3,
            a: [
                [1.3, 2.4, 3.5],
                [2.6, 2.8, 3.3],
                [-1.3, -4.3, -9.7]
            ],
            x: [-2.0, -3.0, 0.0],
            y: [-1.1, 5.0, 0.0],
            incX: 1,
            incY: 1,
            want: [
                [3.5, -7.6, 3.5], 
                [5.9, -12.2, 3.3], 
                [-1.3, -4.3, -9.7]
            ]
        ),
        dgerTestStruct(
            name: "M lt N inc 1",
            m: 3,
            n: 6,
            a: [
				[1.3, 2.4, 3.5, 4.8, 1.11, -9.0],
				[2.6, 2.8, 3.3, -3.4, 6.2, -8.7],
				[-1.3, -4.3, -9.7, -3.1, 8.9, 8.9]
            ],
            x: [-2.0, -3.0, 0.0],
            y: [-1.1, 5.0, 0.0, 9.0, 19.0, 22.0],
            incX: 1,
            incY: 1,
            want: [
                [3.5, -7.6, 3.5, -13.2, -36.89, -53.0],
                [5.9, -12.2, 3.3, -30.4, -50.8, -74.7], 
                [-1.3, -4.3, -9.7, -3.1, 8.9, 8.9]
            ]
        ),
        dgerTestStruct(
            name: "M gt N inc not 1",
            m: 5,
            n: 3,
            a: [
                [1.3, 2.4, 3.5],
                [2.6, 2.8, 3.3],
                [-1.3, -4.3, -9.7],
                [8.0, 9.0, -10.0],
                [-12.0, -14.0, -6.0]
            ],
            x: [-2.0, -3.0, 0.0, 1.0, 2.0, 6.0, 0.0, 9.0, 7.0],
            y: [-1.1, 5.0, 0.0, 8.0, 7.0, -5.0, 7.0],
            incX: 2,
            incY: 3,
            want: [
                [3.5, -13.6, -10.5], 
                [2.6, 2.8, 3.3], 
                [-3.5, 11.7, 4.3], 
                [8.0, 9.0, -10.0], 
                [-19.700000000000003, 42.0, 43.0]
            ]
        ),
        dgerTestStruct(
            name: "M eq N inc not 1",
            m: 3,
            n: 3,
            a: [
                [1.3, 2.4, 3.5],
                [2.6, 2.8, 3.3],
                [-1.3, -4.3, -9.7]
            ],
            x: [-2.0, -3.0, 0.0, 8.0, 7.0, -9.0, 7.0, -6.0, 12.0, 6.0, 6.0, 6.0, -11.0],
            y: [-1.1, 5.0, 0.0, 0.0, 9.0, 8.0, 6.0],
            incX: 4,
            incY: 3,
            want: [
                [3.5, 2.4, -8.5], 
                [-5.1, 2.8, 45.3], 
                [-14.5, -4.3, 62.3]
            ]
        ),
        dgerTestStruct(
            name: "M lt N inc not 1",
            m: 3,
            n: 6,
            a: [
				[1.3, 2.4, 3.5, 4.8, 1.11, -9.0],
				[2.6, 2.8, 3.3, -3.4, 6.2, -8.7],
				[-1.3, -4.3, -9.7, -3.1, 8.9, 8.9]
            ],
            x: [-2.0, -3.0, 0.0, 0.0, 8.0, 0.0, 9.0, -3.0],
            y: [-1.1, 5.0, 0.0, 9.0, 19.0, 22.0, 11.0, -8.11, -9.22, 9.87, 7.0],
            incX: 3,
            incY: 2,
            want: [
                [3.5, 2.4, -34.5, -17.2, 19.55, -23.0], 
                [2.6, 2.8, 3.3, -3.4, 6.2, -8.7], 
                [-11.2, -4.3, 161.3, 95.9, -74.08, 71.9]
            ]
        ),
        dgerTestStruct(
            name: "Y NaN element",
            m: 1,
            n: 1,
            a: [
                [1.3]
            ],
            x: [1.3],
            y: [Float64.NaN],
            incX: 1,
            incY: 1,
            want: [
                [Float64.NaN]
            ]
        ),
        dgerTestStruct(
            name: "M eq N large inc 1",
            m: 7,
            n: 7,
            a: [
				[-20.5, 17.1, -8.4, -23.8, 3.9, 7.7, 6.25],
				[2.9, -0.29, 25.6, -9.4, 36.5, 9.7, 2.3],
				[4.1, -34.1, 10.3, 4.5, -42.05, 9.4, 4.0],
				[19.2, 9.8, -32.7, 4.1, 4.4, -22.5, -7.8],
				[3.6, -24.5, 21.7, 8.6, -13.82, 38.05, -2.29],
				[39.4, -40.5, 7.9, -2.5, -7.7, 18.1, -25.5],
				[-18.5, 43.2, 2.1, 30.1, 3.02, -31.1, -7.6]
            ],
            x: [6.2, -5.0, 88.68, 43.4, -30.5, -40.2, 19.9],
            y: [1.5, 21.7, -28.7, -11.9, 18.1, 3.1, 21.0],
            incX: 1,
            incY: 1,
            want: [
				[-11.2, 151.64, -186.34, -97.58, 116.12, 26.92, 136.45],
				[-4.6, -108.79, 169.1, 50.1, -54.0, -5.8, -102.7],
				[137.12, 1890.256, -2534.816, -1050.792, 1563.058, 284.308, 1866.28],
				[84.3, 951.58, -1278.28, -512.36, 789.94, 112.04, 903.6],
				[-42.15, -686.35, 897.05, 371.55, -565.87, -56.5, -642.79],
				[-20.9, -912.84, 1161.64, 475.88, -735.32, -106.52, -869.7],
				[11.35, 475.03, -569.03, -206.71, 363.21, 30.59, 410.3]
            ]
        ),
        dgerTestStruct(
            name: "M eq N large inc not 1",
            m: 7,
            n: 7,
            a: [
				[-20.5, 17.1, -8.4, -23.8, 3.9, 7.7, 6.25],
				[2.9, -0.29, 25.6, -9.4, 36.5, 9.7, 2.3],
				[4.1, -34.1, 10.3, 4.5, -42.05, 9.4, 4.0],
				[19.2, 9.8, -32.7, 4.1, 4.4, -22.5, -7.8],
				[3.6, -24.5, 21.7, 8.6, -13.82, 38.05, -2.29],
				[39.4, -40.5, 7.9, -2.5, -7.7, 18.1, -25.5],
				[-18.5, 43.2, 2.1, 30.1, 3.02, -31.1, -7.6]
            ],
            x: [6.2, 100.0, 200.0, -5.0, 300.0, 400.0, 88.68, 100.0, 200.0, 43.4, 300.0, 400.0, -30.5, 100.0, 200.0, -40.2, 300.0, 400.0, 19.9],
            y: [1.5, 100.0, 200.0, 300.0, 21.7, 100.0, 200.0, 300.0, -28.7, 100.0, 200.0, 300.0, -11.9, 100.0, 200.0, 300.0, 18.1, 100.0, 200.0, 300.0, 3.1, 100.0, 200.0, 300.0, 21.0],
            incX: 3,
            incY: 4,
            want: [
				[-11.2, 151.64, -186.34, -97.58, 116.12, 26.92, 136.45],
				[-4.6, -108.79, 169.1, 50.1, -54.0, -5.8, -102.7],
				[137.12, 1890.256, -2534.816, -1050.792, 1563.058, 284.308, 1866.28],
				[84.3, 951.58, -1278.28, -512.36, 789.94, 112.04, 903.6],
				[-42.15, -686.35, 897.05, 371.55, -565.87, -56.5, -642.79],
				[-20.9, -912.84, 1161.64, 475.88, -735.32, -106.52, -869.7],
				[11.35, 475.03, -569.03, -206.71, 363.21, 30.59, 410.3]
            ]
        )
    ]
    for (i in 0..tests.size) {
        let test = tests[i]
        let x = sliceCopy(test.x)
        let y = sliceCopy(test.y)

        let a = sliceOfSliceCopy(test.a)

        var alpha = 1.0
        var aFlat = flatten(a)
        nFloat64.dger(test.m, test.n, alpha, x, test.incX, y, test.incY, aFlat, test.n)
        var ans = unflatten(aFlat, test.m, test.n)
        dgercomp(x, test.x, y, test.y, ans, test.want, test.name+" row maj")

        alpha = 4.0
        aFlat = flatten(a)
        nFloat64.dger(test.m, test.n, alpha, x, test.incX, y, test.incY, aFlat, test.n)
        ans = unflatten(aFlat, test.m, test.n)
        let trueCopy = sliceOfSliceCopy(test.want)
        for (i in 0..trueCopy.size) {
            for (j in 0..trueCopy[i].size) {
                trueCopy[i][j] = alpha * (trueCopy[i][j] - a[i][j]) + a[i][j]
            }
        }
        dgercomp(x, test.x, y, test.y, ans, trueCopy, test.name+" row maj alpha")
    }
    return true
}

func dgercomp(x: Array<Float64>, xCopy: Array<Float64>, y: Array<Float64>, yCopy: Array<Float64>, ans: Array<Array<Float64>>, trueAns: Array<Array<Float64>>, name: String): Unit {
    if (!dSliceEqual(x, xCopy)) {
        throw Exception("case ${name}: x modified during call to dger\n%${x}\n${xCopy}")
    }
    if (!dSliceEqual(y, yCopy)) {
        throw Exception("case ${name}: x modified during call to dger\n%${y}\n${yCopy}")
    }

    for (i in 0..ans.size) {
        if (!dSliceTolEqual(ans[i], trueAns[i])) {
            throw Exception("case ${name}: answer mismatch at ${i}.\nExpected ${trueAns},\nFound ${ans}")
        }
    }
}