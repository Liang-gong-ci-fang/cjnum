package cjnum.blas.blas32

import cjnum.blas.*
import std.math.*
import std.unittest.*
import std.unittest.testmacro.*

func newGeneralFrom(a: GeneralCols): General {
    let t = General(rows: a.rows, cols: a.cols, stride: a.cols, data: Array<Float32>(a.rows * a.cols, item: 0.0))
    t.from(a)
    return t
}

func newGeneralColsFrom(a: General): GeneralCols {
    let t = GeneralCols(rows: a.rows, cols: a.cols, stride: a.rows, data: Array<Float32>(a.rows * a.cols, item: 0.0))
    t.from(a)
    return t
}

let generalTests = [
    General(rows: 2, cols: 3, stride: 3, data: [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]),
    General(rows: 3, cols: 2, stride: 2, data: [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]),
    General(rows: 3, cols: 3, stride: 3, data: [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]),
    General(rows: 2, cols: 3, stride: 5, data: [1.0, 2.0, 3.0, 0.0, 0.0, 4.0, 5.0, 6.0, 0.0, 0.0]),
    General(rows: 3, cols: 2, stride: 5, data: [1.0, 2.0, 0.0, 0.0, 0.0, 3.0, 4.0, 0.0, 0.0, 0.0, 5.0, 6.0, 0.0, 0.0,
        0.0]),
    General(rows: 3, cols: 3, stride: 5, data: [1.0, 2.0, 3.0, 0.0, 0.0, 4.0, 5.0, 6.0, 0.0, 0.0, 7.0, 8.0, 9.0, 0.0,
        0.0])
]

@Test
public func testConvertGeneral() {
    for (test in generalTests) {
        let colmajor = newGeneralColsFrom(test)
        @Expect(sameGeneral(colmajor, test))
        let rowmajor = newGeneralFrom(colmajor)
        @Expect(sameGeneral(rowmajor, test))
    }
}

extend General <: general {
    public func dims(): (Int64, Int64) {
        return (this.rows, this.cols)
    }

    public func at(i: Int64, j: Int64): Float32 {
        return this.data[i * this.stride + j]
    }
}

extend GeneralCols <: general {
    public func dims(): (Int64, Int64) {
        return (this.rows, this.cols)
    }

    public func at(i: Int64, j: Int64): Float32 {
        return this.data[i + j * this.stride]
    }
}

sealed interface general {
    func dims(): (Int64, Int64)
    func at(i: Int64, j: Int64): Float32
}

func sameGeneral(a: general, b: general): Bool {
    let (ar, ac) = a.dims()
    let (br, bc) = b.dims()
    if (ar != br || ac != bc) {
        return false
    }
    for (i in 0..ar) {
        for (j in 0..ac) {
            if (a.at(i, j) != b.at(i, j) || a.at(i, j).isNaN() != b.at(i, j).isNaN()) {
                return false
            }
        }
    }
    return true
}
