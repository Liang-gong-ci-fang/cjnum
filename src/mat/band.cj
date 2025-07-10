package cjnum.mat

import cjnum.blas.blas64.*
import std.math.*

/*
因为以下注释代码本质上是为了在编译过程中检查是BandDense类否实现了以下接口，
仓颉编译过程中会对未实现的接口方法进行提示，
所以这里没有实现接口判断
 */
// TODO var (
// 	bandDense *BandDense
// 	_         Matrix      = bandDense
// 	_         allMatrix   = bandDense
// 	_         denseMatrix = bandDense
// 	_         Banded      = bandDense
// 	_         RawBander   = bandDense
     
// 	_ NonZeroDoer    = bandDense
// 	_ RowNonZeroDoer = bandDense
// 	_ ColNonZeroDoer = bandDense
// )

// BandDense represents a band matrix in dense storage format.
public class BandDense <: Matrix & DiagDenseViewer {
    let mat: Band
    public init(mat: Band) {
        this.mat = mat
    }

    // Dims returns the number of rows and columns in the matrix.
    public func dims(): (Int64, Int64) {
        return (this.mat.rows, this.mat.cols)
    }

    // Bandwidth returns the upper and lower bandwidths of the matrix.
    public func bandwidth(): (Int64, Int64) {
        return (this.mat.kl, this.mat.ku)
    }

    public func t(): Matrix {
        return Transpose(this)
    }

    // At returns the element at row i, column j.
    @When[feature != "bounds"]
    public func at(i: Int64, j: Int64): Float64 {
        if (i < 0 || UInt64(i) >= UInt64(this.mat.rows)) {
            throw Exception(errRowAccess)
        }
        if (j < 0 || UInt64(j) >= UInt64(this.mat.cols)) {
            throw Exception(errColAccess)
        }
        return this._at(i, j)
    }

    @When[feature != "bounds"]
    private func _at(i: Int64, j: Int64): Float64 {
        let pj = j + this.mat.kl - i
        if (pj < 0 || this.mat.kl + this.mat.ku + 1 <= pj) {
            return 0.0
        }
        return this.mat.data[i * this.mat.stride + pj]
    }

    // At returns the element at row i, column j.
    @When[feature == "bounds"]
    public func at(i: Int64, j: Int64): Float64 {
        return this._at(i, j)
    }

    @When[feature == "bounds"]
    private func _at(i: Int64, j: Int64): Float64 {
        if (i < 0 || UInt64(i) >= UInt64(this.mat.rows)) {
            throw Exception(errRowAccess)
        }
        if (j < 0 || UInt64(j) >= UInt64(this.mat.cols)) {
            throw Exception(errColAccess)
        }
        let pj = j + this.mat.kl - i
        if (pj < 0 || this.mat.kl + this.mat.ku + 1 <= pj) {
            return 0.0
        }
        return this.mat.data[i * this.mat.stride + pj]
    }

    // SetBand sets the element at row i, column j to the value v.
    // It panics if the location is outside the appropriate region of the matrix.
    @When[feature != "bounds"]
    public func setBand(i: Int64, j: Int64, v: Float64): Unit {
        if (i < 0 || UInt64(i) >= UInt64(this.mat.rows)) {
            throw Exception(errRowAccess)
        }
        if (j < 0 || UInt64(j) >= UInt64(this.mat.cols)) {
            throw Exception(errColAccess)
        }

        let pj = j + this.mat.kl - i
        if (pj < 0 || this.mat.kl + this.mat.ku + 1 <= pj) {
            throw Exception(errBandSet)
        }

        this.set(i, j, v)
    }

    @When[feature != "bounds"]
    private func set(i: Int64, j: Int64, v: Float64): Unit {
        let pj = j + this.mat.kl - i

        this.mat.data[i * this.mat.stride + pj] = v
    }

    // SetBand sets the element at row i, column j to the value v.
    // It panics if the location is outside the appropriate region of the matrix.
    @When[feature == "bounds"]
    public func setBand(i: Int64, j: Int64, v: Float64) {
        this.set(i, j, v)
    }

    @When[feature == "bounds"]
    private func set(i: Int64, j: Int64, v: Float64) {
        if (i < 0 || UInt64(i) >= UInt64(this.mat.rows)) {
            throw Exception(errRowAccess)
        }
        if (j < 0 || UInt64(j) >= UInt64(this.mat.cols)) {
            throw Exception(errColAccess)
        }
        let pj = j + this.mat.kl - i
        if (pj < 0 || this.mat.kl + this.mat.ku + 1 <= pj) {
            throw Exception(errBandSet)
        }
        this.mat.data[i * this.mat.stride + pj] = v
    }

    // Zero sets all of the matrix elements to zero.
    public func zero() {
        let m = this.mat.rows
        let kL = this.mat.kl
        let nCol = this.mat.ku + 1 + kL
        for (i in 0..m) {
            let l = max(0, kL - i)
            let u = min(nCol, m + kL - i)
            _zero(this.mat.data[i * this.mat.stride + l..i * this.mat.stride + u])
        }
    }

    // DiagView returns the diagonal as a matrix backed by the original data.
    public func diagView(): Diagonal {
        let n = min(this.mat.rows, this.mat.cols)
        return DiagDense(
            Vector(
                n: n,
                inc: this.mat.stride,
                data: this.mat.data[this.mat.kl..(n - 1) * this.mat.stride + this.mat.kl + 1],
            ),
        )
    }
}

let bendDense = BandDense(
    Band(
        rows: 6,
        cols: 6,
        kl: 1,
        ku: 2,
        stride: 4,
        data: [-1.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0,
            -1.0, 19.0, 20.0, -1.0, -1.0]
    )
)
let judge = bendDense is Matrix

// A RawBander can return a blas64.Band representation of the receiver.
// Changes to the blas64.Band.Data slice will be reflected in the original
// matrix, changes to the Rows, Cols, KL, KU and Stride fields will not.
public interface RawBander {
    func rawBand(): Band
}

// NewBandDense creates a new Band matrix with r rows and c columns. If data == nil,
// a new slice is allocated for the backing slice. If len(data) == min(r, c+kl)*(kl+ku+1),
// data is used as the backing slice, and changes to the elements of the returned
// BandDense will be reflected in data. If neither of these is true, NewBandDense
// will panic. kl must be at least zero and less r, and ku must be at least zero and
// less than c, otherwise NewBandDense will panic.
// NewBandDense will panic if either r or c is zero.
//
// The data must be arranged in row-major order constructed by removing the zeros
// from the rows outside the band and aligning the diagonals. For example, the matrix
//
//	1  2  3  0  0  0
//	4  5  6  7  0  0
//	0  8  9 10 11  0
//	0  0 12 13 14 15
//	0  0  0 16 17 18
//	0  0  0  0 19 20
//
// becomes (* entries are never accessed)
//   - 1  2  3
//     4  5  6  7
//     8  9 10 11
//     12 13 14 15
//     16 17 18  *
//     19 20  *  *
//
// which is passed to NewBandDense as []float64{*, 1, 2, 3, 4, ...} with kl=1 and ku=2.
// Only the values in the band portion of the matrix are used.
public func newBandDense(r: Int64, c: Int64, kl: Int64, ku: Int64, data_: Array<Float64>): BandDense {
    if (r <= 0 || c <= 0 || kl < 0 || ku < 0) {
        if (r == 0 || c == 0) {
            throw Exception(errZeroLength)
        }
        throw Exception(errNegativeDimension)
    }
    if (kl + 1 > r || ku + 1 > c) {
        throw Exception(errBandwidth)
    }
    let bc = kl + ku + 1
    var data = data_
    if (data.size > 0 && data.size != min(r, c + kl) * bc) {
        throw Exception(errShape)
    }
    if (data.size == 0) {
        data = Array<Float64>(min(r, c + kl) * bc, item: 0.0)
    }
    return BandDense(
        Band(
            rows: r,
            cols: c,
            kl: kl,
            ku: ku,
            stride: bc,
            data: data,
        ),
    )
}

// NewDiagonalRect is a convenience function that returns a diagonal matrix represented by a
// BandDense. The length of data must be min(r, c) otherwise NewDiagonalRect will panic.
public func newDiagonalRect(r: Int64, c: Int64, data: Array<Float64>): BandDense {
    return newBandDense(r, c, 0, 0, data)
}
