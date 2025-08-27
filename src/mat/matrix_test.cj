package cjnum.mat

import std.core.*
import std.collection.*

public func mPanics(fn: () -> Unit): (Bool, String) {
    try {
        fn()
    } catch (e: Exception) {
        return (true, e.message)
    }
    return (false, "")
}

public func flatten(f: Array<Array<Float64>>): (Int64, Int64, Array<Float64>) {
    let r = f.size
    if (r == 0) {
        throw Exception("bad test: no row")
    }
    let c = f[0].size
    let d = ArrayList<Float64>(r * c)
    for (row in f) {
        if (row.size != c) {
            throw Exception("bad test: ragged input")
        }
        d.add(all: row)
    }
    // println("d.size = ${d.size}")
    return (r, c, d.toArray())
}

//TODO @Test
// public func testCol(): Unit {
//     let colTests: Array<Array<Array<Float64>>> = [
//         [
//             [1.0, 2.0, 3.0],
//             [4.0, 5.0, 6.0],
//             [7.0, 8.0, 9.0]
//         ],
//         [
//             [1.0, 2.0, 3.0],
//             [4.0, 5.0, 6.0],
//             [7.0, 8.0, 9.0],
//             [10.0, 11.0, 12.0]
//         ],
//         [
//             [1.0, 2.0, 3.0, 4.0],
//             [5.0, 6.0, 7.0, 8.0],
//             [9.0, 10.0, 11.0, 12.0]
//         ]
//     ]
//     for (id in 0..colTests.size) {
//         let af = colTests[id]
//         let (af_1, af_2, af_3) = flatten(af)
//         let a = newDense(af_1, af_2, af_3)
//         var col = Array<Float64>(a.mat.rows, repeat: 0.0)
//         for (j in 0..af[0].size) {
//             for (i in 0..col.size) {
//                 col[i] = Float64(i * a.mat.cols + j + 1)
//             }

//             var got = mCol(Array<Float64>(0, repeat: 0.0), j, a)
//             @Expect(got==col)

//             got = Array<Float64>(a.mat.rows, repeat: 0.0)
//             mCol(got, j, a)
//             @Expect(got==col)
//         }
//     }

//     func denseComparison(a: Dense): Any {
//         let (r, c) = a.dims()
//         var ans = ArrayList<Array<Float64>>(c)
//         for (j in 0..ans.size) {
//             ans[j] = Array<Float64>(r, repeat: 0.0)
//             for (i in 0..ans[j].size) {
//                 ans[j][i] = a.at(i, j)
//             }
//         }
//         return ans
//     }

//     func f(a: Matrix): Any {
//         let (_, c) = a.dims()
//         var ans = ArrayList<Array<Float64>>(c)
//         for (j in 0..ans.size) {
//             ans[j] = mCol(Array<Float64>(0, repeat: 0.0), j, a)
//         }
//         return ans
//     }
//     testOneInputFunc(t, "Col", f, denseComparison, sameAnswerF64SliceOfSlice, isAnyType, isAnySize)

//     // f = func(a Matrix) interface{} {
//     // 	r, c := a.Dims()
//     // 	ans := make([][]float64, c)
//     // 	for j := range ans {
//     // 		ans[j] = make([]float64, r)
//     // 		Col(ans[j], j, a)
//     // 	}
//     // 	return ans
//     // }
//     // testOneInputFunc(t, "Col", f, denseComparison, sameAnswerF64SliceOfSlice, isAnyType, isAnySize)
// }
