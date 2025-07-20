package cjnum.blas.blas64

import cjnum.blas.*
import std.math.*
import std.unittest.*
import std.unittest.testmacro.*

func newBandFrom(a: BandCols): Band {
    let t = Band(rows: a.rows, cols: a.cols, kl: a.kl, ku: a.ku, stride: a.kl + a.ku + 1,
                 data: Array<Float64>(a.rows * (a.kl + a.ku + 1), repeat: 0.0))
    for (i in 0..t.data.size) {
        t.data[i] = Float64.NaN
    }
    t.from(a)
    return t
}

extend Band <: band {
    public func dims(): (Int64, Int64) {
        return (this.rows, this.cols)
    }

    public func at(i: Int64, j: Int64): Float64 {
        let pj = j + this.kl - i
        if (pj < 0 || this.kl + this.ku + 1 <= pj) {
            return 0.0
        }
        return this.data[i * this.stride + pj]
    }

    public func bandwidth(): (Int64, Int64) {
        return (this.kl, this.ku)
    }
}

func newBandColsFrom(a: Band): BandCols {
let t = BandCols(rows: a.rows, cols: a.cols, kl: a.kl, ku: a.ku, stride: a.kl + a.ku + 1,
                 data: Array<Float64>(a.cols * (a.kl + a.ku + 1), repeat: 0.0))
    for (i in 0..t.data.size) {
        t.data[i] = Float64.NaN
    }
    t.from(a)
    return t
}

extend BandCols <: band {
    public func dims(): (Int64, Int64) {
        return (this.rows, this.cols)
    }

    public func at(i: Int64, j: Int64): Float64 {
        let pj = i + this.ku - j
        if (pj < 0 || this.kl + this.ku + 1 <= pj) {
            return 0.0
        }
        return this.data[j * this.stride + pj]
    }

    public func bandwidth(): (Int64, Int64) {
        return (this.kl, this.ku)
    }
}

sealed interface band {
    func dims(): (Int64, Int64)
    func at(i: Int64, j: Int64): Float64
    func bandwidth(): (Int64, Int64)
}

func sameBand(a: band, b: band): Bool {
    let (ar, ac) = a.dims()
    let (br, bc) = b.dims()
    if (ar != br || ac != bc) {
        return false
    }
    let (akl, aku) = a.bandwidth()
    let (bkl, bku) = b.bandwidth()
    if (akl != bkl || aku != bku) {
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

let bandTests = [
    Band(rows: 3, cols: 4, kl: 0, ku: 0, stride: 1, data: [1.0, 2.0, 3.0]),
    Band(rows: 3, cols: 3, kl: 0, ku: 0, stride: 1, data: [1.0, 2.0, 3.0]),
    Band(rows: 4, cols: 3, kl: 0, ku: 0, stride: 1, data: [1.0, 2.0, 3.0]),
    Band(rows: 4, cols: 3, kl: 0, ku: 1, stride: 2, data: [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]),
    Band(rows: 3, cols: 4, kl: 0, ku: 1, stride: 2, data: [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]),
    Band(rows: 3, cols: 4, kl: 1, ku: 1, stride: 3, data: [-1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]),
    Band(rows: 4, cols: 3, kl: 1, ku: 1, stride: 3, data: [-1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, -2.0, 9.0, -3.0, -4.0]),
    Band(rows: 3, cols: 4, kl: 2, ku: 1, stride: 4, data: [-2.0, -1.0, 3.0, 4.0, -3.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0]),
    Band(rows: 4, cols: 3, kl: 2, ku: 1, stride: 4, data: [-2.0, -1.0, 2.0, 3.0, -3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, -4.0, 10.0, 11.0,
    -5.0, -6.0]),
    Band(rows: 3, cols: 4, kl: 0, ku: 0, stride: 5, data: [1.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 0.0]),
    Band(rows: 3, cols: 3, kl: 0, ku: 0, stride: 5, data: [1.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 0.0]),
    Band(rows: 4, cols: 3, kl: 0, ku: 0, stride: 5, data: [1.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 0.0]),
    Band(rows: 4, cols: 3, kl: 0, ku: 1, stride: 5, data: [1.0, 2.0, 0.0, 0.0, 0.0, 3.0, 4.0, 0.0, 0.0, 0.0, 5.0, 6.0, 0.0, 0.0, 0.0]),
    Band(rows: 3, cols: 4, kl: 0, ku: 1, stride: 5, data: [1.0, 2.0, 0.0, 0.0, 0.0, 3.0, 4.0, 0.0, 0.0, 0.0, 5.0, 6.0, 0.0, 0.0, 0.0]),
    Band(rows: 3, cols: 4, kl: 1, ku: 1, stride: 5, data: [-1.0, 2.0, 3.0, 0.0, 0.0, 4.0, 5.0, 6.0, 0.0, 0.0, 7.0, 8.0, 9.0, 0.0, 0.0]),
    Band(rows: 4, cols: 3, kl: 1, ku: 1, stride: 5, data: [-1.0, 2.0, 3.0, 0.0, 0.0, 4.0, 5.0, 6.0, 0.0, 0.0, 7.0, 8.0, -2.0, 0.0, 0.0, 
    9.0, -3.0, -4.0, 0.0, 0.0]),
    Band(rows: 3, cols: 4, kl: 2, ku: 1, stride: 5, data: [-2.0, -1.0, 3.0, 4.0, 0.0, -3.0, 5.0, 6.0, 7.0, 0.0, 8.0, 9.0, 10.0, 11.0, 0.0]),
    Band(rows: 4, cols: 3, kl: 2, ku: 1, stride: 5, data: [-2.0, 1.0, 2.0, 3.0, 0.0, -3.0, 4.0, 5.0, 6.0, 0.0, 7.0, 8.0, 9.0, -4.0, 0.0, 
    10.0, 11.0, -5.0, -6.0, 0.0])
]

@Test
public func testConvertBand() {
    for (test in bandTests) {
        let colmajor = newBandColsFrom(test)
        @Expect(sameBand(colmajor, test))
        let rowmajor = newBandFrom(colmajor)
        @Expect(sameBand(rowmajor, test))
    }
}
