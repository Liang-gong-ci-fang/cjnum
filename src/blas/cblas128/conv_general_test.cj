package cjnum.blas.cblas128

import cjnum.complex.*
import std.math.*
import std.unittest.*
import std.unittest.testmacro.*

func newGeneralFrom(a: GeneralCols): General {
    let t = General(rows: a.rows, cols: a.cols, stride: a.cols, data: Array<Complex128>(a.rows * a.cols, repeat: Complex128(0.0)))
    t.from(a)
    return t
}

func newGeneralColsFrom(a: General): GeneralCols {
    let t = GeneralCols(rows: a.rows, cols: a.cols, stride: a.rows, data: Array<Complex128>(a.rows * a.cols, repeat: Complex128(0.0)))
    t.from(a)
    return t
}

let generalTests = [
    General(
        rows: 2, 
        cols: 3, 
        stride: 3,
        data: [
            Complex128(1.0), Complex128(2.0), Complex128(3.0),
            Complex128(4.0), Complex128(5.0), Complex128(6.0)
            ]
        ),
    General(
        rows: 3, 
        cols: 2, 
        stride: 2, 
        data: [
            Complex128(1.0), Complex128(2.0), 
            Complex128(3.0), Complex128(4.0), 
            Complex128(5.0), Complex128(6.0)
            ]
        ),
    General(
        rows: 3, 
        cols: 3,
        stride: 3, 
        data: [
            Complex128(1.0), Complex128(2.0), Complex128(3.0), 
            Complex128(4.0), Complex128(5.0), Complex128(6.0), 
            Complex128(7.0), Complex128(8.0), Complex128(9.0)
            ]
        ),
    General(
        rows: 2, 
        cols: 3, 
        stride: 5, 
        data: [
            Complex128(1.0), Complex128(2.0), Complex128(3.0), Complex128(0.0), Complex128(0.0), 
            Complex128(4.0), Complex128(5.0), Complex128(6.0), Complex128(0.0), Complex128(0.0)
            ]
        ),
    General(
        rows: 3, 
        cols: 2, 
        stride: 5, 
        data: [
            Complex128(1.0), Complex128(2.0), Complex128(0.0), Complex128(0.0), Complex128(0.0), 
            Complex128(3.0), Complex128(4.0), Complex128(0.0), Complex128(0.0), Complex128(0.0), 
            Complex128(5.0), Complex128(6.0), Complex128(0.0), Complex128(0.0), Complex128(0.0)
            ]
        ),
    General(
        rows: 3, 
        cols: 3, 
        stride: 5, 
        data: [
            Complex128(1.0), Complex128(2.0), Complex128(3.0), Complex128(0.0), Complex128(0.0), 
            Complex128(4.0), Complex128(5.0), Complex128(6.0), Complex128(0.0), Complex128(0.0), 
            Complex128(7.0), Complex128(8.0), Complex128(9.0), Complex128(0.0), Complex128(0.0)
            ]
        )
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

    public func at(i: Int64, j: Int64): Complex128 {
        return this.data[i * this.stride + j]
    }
}

extend GeneralCols <: general {
    public func dims(): (Int64, Int64) {
        return (this.rows, this.cols)
    }

    public func at(i: Int64, j: Int64): Complex128 {
        return this.data[i + j * this.stride]
    }
}

sealed interface general {
    func dims(): (Int64, Int64)
    func at(i: Int64, j: Int64): Complex128
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
