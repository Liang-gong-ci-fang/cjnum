package cjnum.mat

import cjnum.blas.*
import cjnum.blas.blas64.*

/*
因为以下注释代码本质上是为了在编译过程中检查是SymBandDense类否实现了以下接口，
仓颉编译过程中会对未实现的接口方法进行提示，
所以这里没有实现接口判断
*/
// TODO var (
// 	symBandDense *SymBandDense
// 	_            Matrix           = symBandDense
// 	_            allMatrix        = symBandDense
// 	_            denseMatrix      = symBandDense
// 	_            Symmetric        = symBandDense
// 	_            Banded           = symBandDense
// 	_            SymBanded        = symBandDense
// 	_            RawSymBander     = symBandDense
// 	_            MutableSymBanded = symBandDense

// 	_ NonZeroDoer    = symBandDense
// 	_ RowNonZeroDoer = symBandDense
// 	_ ColNonZeroDoer = symBandDense
// )


// SymBandDense represents a symmetric band matrix in dense storage format.
public class SymBandDense <: Matrix {
    let mat: SymmetricBand
    public init(mat: SymmetricBand) {
        this.mat = mat
    }

    @When[feature != "bounds"]
    public func at(i: Int64, j: Int64): Float64 {
        if (i < 0 || UInt64(i) >= UInt64(s.mat.n)) {
            throw Exception(errRowAccess)
        }
        if (j < 0 || UInt64(j) >= UInt64(s.mat.N)) {
            throw Exception(errColAccess)
        }
        return s._at(i, j)
    }

    @When[feature != "bounds"]
    private func _at(i: Int64, j: Int64): Float64 {
        if (i > j) {
            let i_ = j
            let j_ = i
        }
        let pj = j - i
        if (s.mat.K + 1 <= pj) {
            return 0
        }
        return s.mat.Data[i * s.mat.Stride + pj]
    }

    // At returns the element at row i, column j.
    @When[feature == "bounds"]
    public func at(i: Int64, j: Int64): Float64 {
        return this._at(i, j)
    }

    @When[feature == "bounds"]
    private func _at(i: Int64, j: Int64): Float64 {
        if (i < 0 || UInt64(i) >= UInt64(this.mat.n)) {
            throw Exception(errRowAccess)
        }
        if (j < 0 || UInt64(j) >= UInt64(this.mat.n)) {
            throw Exception(errColAccess)
        }
        if (i > j) {
            let i_ = j
            let j_ = i
            let pj = j_ - i_
            if (this.mat.k + 1 <= pj) {
                return 0.0
            }
            return this.mat.data[i_ * this.mat.stride + pj]
        }
        let pj = j - i
        if (this.mat.k + 1 <= pj) {
            return 0.0
        }
        return this.mat.data[i * this.mat.stride + pj]
    }

    // Dims returns the number of rows and columns in the matrix.
    public func dims(): (Int64, Int64) {
        return (this.mat.n, this.mat.n)
    }

    // T implements the Matrix interface. Symmetric matrices, by definition, are
    // equal to their transpose, and this is a no-op.
    public func t(): Matrix {
        return this
    }
}

// A RawSymBander can return a blas64.SymmetricBand representation of the receiver.
// Changes to the blas64.SymmetricBand.Data slice will be reflected in the original
// matrix, changes to the N, K, Stride and Uplo fields will not.
public interface RawSymBander {
    func rawSymBand(): SymmetricBand
}

// NewSymBandDense creates a new SymBand matrix with n rows and columns. If data == nil,
// a new slice is allocated for the backing slice. If len(data) == n*(k+1),
// data is used as the backing slice, and changes to the elements of the returned
// SymBandDense will be reflected in data. If neither of these is true, NewSymBandDense
// will panic. k must be at least zero and less than n, otherwise NewSymBandDense will panic.
//
// The data must be arranged in row-major order constructed by removing the zeros
// from the rows outside the band and aligning the diagonals. SymBandDense matrices
// are stored in the upper triangle. For example, the matrix
//
//	1  2  3  0  0  0
//	2  4  5  6  0  0
//	3  5  7  8  9  0
//	0  6  8 10 11 12
//	0  0  9 11 13 14
//	0  0  0 12 14 15
//
// becomes (* entries are never accessed)
//
//	 1  2  3
//	 4  5  6
//	 7  8  9
//	10 11 12
//	13 14  *
//	15  *  *
//
// which is passed to NewSymBandDense as []float64{1, 2, ..., 15, *, *, *} with k=2.
// Only the values in the band portion of the matrix are used.
public func newSymBandDense(n: Int64, k: Int64, data: Array<Float64>): SymBandDense {
    if (n <= 0 || k < 0) {
        if (n == 0) {
            throw Exception(errZeroLength)
        }
        throw Exception("mat: negative dimension")
    }
    if (k + 1 > n) {
        throw Exception("mat: band out of range")
    }
    let bc = k + 1
    if (data.size != 0 && data.size != n * bc) {
        throw Exception(errShape)
    }
    if (data.size == 0) {
        let data_ = Array<Float64>(n * bc, repeat: 0.0)
        return SymBandDense(
            SymmetricBand(
                n: n,
                k: k,
                stride: bc,
                uplo: Upper,
                data: data_,
            )
        )
    }
    return SymBandDense(
        SymmetricBand(
            n: n,
            k: k,
            stride: bc,
            uplo: Upper,
            data: data,
        )
    )
}
