package cjnum.mat

import cjnum.blas.blas64.*

/*
因为以下注释代码本质上是为了在编译过程中检查是VecDense类否实现了以下接口，
仓颉编译过程中会对未实现的接口方法进行提示，
所以这里没有实现接口判断
 */
// var (
// 	vector *VecDense

// 	_ Matrix        = vector
// 	_ allMatrix     = vector
// 	_ Vector        = vector
// 	_ Reseter       = vector
// 	_ MutableVector = vector
// )

// VecDense represents a column vector.
public class VecDense <: Matrix {
    let mat: Vector
    // A BLAS vector can have a negative increment, but allowing this
    // in the mat type complicates a lot of code, and doesn't gain anything.
    // VecDense must have positive increment in this package.

    public init(mat: Vector) {
        this.mat = mat
    }

    // At returns the element at row i.
    // It panics if i is out of bounds or if j is not zero.
    @When[feature != "bounds"]
    public func at(i: Int64, j: Int64): Float64 {
        if (i < 0 || UInt64(i) >= UInt64(this.mat.N)) {
            throw Exception(ErrRowAccess)
        }
        if (j != 0) {
            throw Exception(ErrColAccess)
        }
        return this._at(i)
    }

    @When[feature != "bounds"]
    public func _at(i: Int64): Float64 {
        return this.mat.data[i * this.mat.inc]
    }

    // At returns the element at row i.
    // It panics if i is out of bounds or if j is not zero.
    @When[feature == "bounds"]
    public func at(i: Int64, j: Int64): Float64 {
        if (j != 0) {
            throw Exception(errColAccess)
        }
        return this._at(i)
    }

    @When[feature == "bounds"]
    public func _at(i: Int64): Float64 {
        if (i < 0 || UInt64(i) >= UInt64(this.mat.n)) {
            throw Exception(errRowAccess)
        }
        return this.mat.data[i * this.mat.inc]
    }

    // Dims returns the number of rows and columns in the matrix. Columns is always 1
    // for a non-Reset vector.
    public func dims(): (Int64, Int64) {
        if (this.isEmpty()) {
            return (0, 0)
        }
        return (this.mat.n, 1)
    }

    // T performs an implicit transpose by returning the receiver inside a Transpose.
    public func t(): Matrix {
        return Transpose(this)
    }

    // IsEmpty returns whether the receiver is empty. Empty matrices can be the
    // receiver for size-restricted operations. The receiver can be emptied using
    // Reset.
    public func isEmpty(): Bool {
        // It must be the case that v.Dims() returns
        // zeros in this case. See comment in Reset().
        return this.mat.inc == 0
    }
}

// NewVecDense creates a new VecDense of length n. If data == nil,
// a new slice is allocated for the backing slice. If len(data) == n, data is
// used as the backing slice, and changes to the elements of the returned VecDense
// will be reflected in data. If neither of these is true, NewVecDense will panic.
// NewVecDense will panic if n is zero.
public func newVecDense(n: Int64, data: Array<Float64>): VecDense {
    if (n <= 0) {
        if (n == 0) {
            throw Exception(errZeroLength)
        }
        throw Exception("mat: negative dimension")
    }
    if (data.size != n && data.size != 0) {
        throw Exception(errShape)
    }
    if (data.size == 0) {
        let data_ = Array<Float64>(n, item: 0.0)
        return VecDense(
            Vector(
                n: n,
                inc: 1,
                data: data_,
            ),
        )
    }
    return VecDense(
        Vector(
            n: n,
            inc: 1,
            data: data,
        ),
    )
}
