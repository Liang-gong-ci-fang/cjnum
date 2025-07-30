package cjnum.mat

import cjnum.blas.*
import cjnum.blas.blas64.*

/*
因为以下注释代码本质上是为了在编译过程中检查是SymDense类否实现了以下接口，
仓颉编译过程中会对未实现的接口方法进行提示，
所以这里没有实现接口判断
*/
// TODO var (
// 	symDense *SymDense

// 	_ Matrix           = symDense
// 	_ allMatrix        = symDense
// 	_ denseMatrix      = symDense
// 	_ Symmetric        = symDense
// 	_ RawSymmetricer   = symDense
// 	_ MutableSymmetric = symDense
// )

const badSymTriangle = "mat: blas64.Symmetric not upper"

// SymDense is a symmetric matrix that uses dense storage. SymDense
// matrices are stored in the upper triangle.
public class SymDense <: Matrix {
    let mat: Symmetric
    let cap: Int64

    public init(mat: Symmetric, cap: Int64) {
        this.mat = mat
        this.cap = cap
    }

    public init(mat: Symmetric) {
        if (mat.uplo != Upper) {
            throw Exception(badSymTriangle)
        }
        this.mat = mat
        this.cap = mat.n
    }

    // At returns the element at row i and column j.
    @When[feature != "bounds"]
    public func at(i: Int64, j: Int64): Float64 {
        if (i < 0 || UInt64(i) >= UInt64(s.mat.N)) {
            throw Exception(errRowAccess)
        }
        if (j < 0 || UInt64(j) >= UInt64(s.mat.N)) {
            throw Exception(errColAccess)
        }
        return s._at(i, j)
    }

    @When[feature != "bounds"]
    public func _at(i: Int64, j: Int64): Float64 {
        if (i > j) {
            let i_ = j
            let j_ = i
            return s.mat.Data[i_ * s.mat.Stride + j_]
        }
        return s.mat.Data[i * s.mat.Stride + j]
    }

    // At returns the element at row i and column j.
    @When[feature == "bounds"]
    public func at(i: Int64, j: Int64): Float64 {
        return this._at(i, j)
    }

    @When[feature == "bounds"]
    public func _at(i: Int64, j: Int64): Float64 {
        if (i < 0 || UInt64(i) >= UInt64(this.mat.n)) {
            throw Exception(errRowAccess)
        }
        if (j < 0 || UInt64(j) >= UInt64(this.mat.n)) {
            throw Exception(errColAccess)
        }
        if (i > j) {
            let i_ = j
            let j_ = i
            return this.mat.data[i_ * this.mat.stride + j_]
        }
        return this.mat.data[i * this.mat.stride + j]
    }

    // Dims returns the number of rows and columns in the matrix.
    public func dims(): (Int64, Int64) {
        return (this.mat.n, this.mat.n)
    }

    // T returns the receiver, the transpose of a symmetric matrix.
    public func t(): Matrix {
        return this
    }
}

// A RawSymmetricer can return a view of itself as a BLAS Symmetric matrix.
public interface RawSymmetricer {
    func rawSymmetric(): Symmetric
}

// NewSymDense creates a new Symmetric matrix with n rows and columns. If data == nil,
// a new slice is allocated for the backing slice. If len(data) == n*n, data is
// used as the backing slice, and changes to the elements of the returned SymDense
// will be reflected in data. If neither of these is true, NewSymDense will panic.
// NewSymDense will panic if n is zero.
//
// The data must be arranged in row-major order, i.e. the (i*c + j)-th
// element in the data slice is the {i, j}-th element in the matrix.
// Only the values in the upper triangular portion of the matrix are used.
public func newSymDense(n: Int64, data: Array<Float64>): SymDense {
    if (n <= 0) {
        if (n == 0) {
            throw Exception(errZeroLength)
        }
        throw Exception("mat: negative dimension")
    }
    if (data.size != 0 && n * n != data.size) {
        throw Exception(errShape)
    }
    if (data.size == 0) {
        let data_ = Array<Float64>(n * n, repeat: 0.0)
        return SymDense(
            Symmetric(
                n: n,
                stride: n,
                data: data_,
                uplo: Upper,
            ),
            n,
        )
    }
    return SymDense(
        Symmetric(
            n: n,
            stride: n,
            data: data,
            uplo: Upper,
        ),
        n,
    )
}
