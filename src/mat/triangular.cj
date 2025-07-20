package cjnum.mat

import cjnum.blas.blas64.*
import cjnum.blas.*

/*
因为以下注释代码本质上是为了在编译过程中检查是TriDense类否实现了以下接口，
仓颉编译过程中会对未实现的接口方法进行提示，
所以这里没有实现接口判断
*/
// TODO var (
// 	triDense *TriDense
// 	_        Matrix            = triDense
// 	_        allMatrix         = triDense
// 	_        denseMatrix       = triDense
// 	_        Triangular        = triDense
// 	_        RawTriangular     = triDense
// 	_        MutableTriangular = triDense

// 	_ NonZeroDoer    = triDense
// 	_ RowNonZeroDoer = triDense
// 	_ ColNonZeroDoer = triDense
// )

// TriDense represents an upper or lower triangular matrix in dense storage
// format.
public class TriDense <: Matrix & RawTriangular {
    let mat: Triangular
    let cap: Int64
    public init(mat: Triangular, cap: Int64) {
        this.mat = mat
        this.cap = cap
    }

    // At returns the element at row i, column j.
    @When[feature != "bounds"]
    public func at(i: Int64, j: Int64): Float64 {
        if (i == Int64.NaN || i < 0 || UInt64(i) >= UInt64(this.mat.n)) {
            throw Exception(errRowAccess)
        }
        if (j == Int64.NaN || j < 0 || UInt64(j) >= UInt64(this.mat.n)) {
            throw Exception(errColAccess)
        }
        return this._at(i, j)
    }

    @When[feature != "bounds"]
    public func _at(i: Int64, j: Int64): Float64 {
        let isUpper = this.triKind()
        if ((isUpper && i > j) || (!isUpper && i < j)) {
            return 0.0
        }
        return this.mat.Data[i * this.mat.Stride + j]
    }

    // At returns the element at row i, column j.
    @When[feature == "bounds"]
    public func at(i: Int64, j: Int64): Float64 {
        return this.at(i, j)
    }

    @When[feature == "bounds"]
    public func _at(i: Int64, j: Int64): Float64 {
        if (i < 0 || UInt64(i) >= UInt64(this.mat.n)) {
            throw Exception(errRowAccess)
        }
        if (j < 0 || UInt64(j) >= UInt64(this.mat.n)) {
            throw Exception(errColAccess)
        }
        let isUpper = this.isUpper()
        if ((isUpper && i > j) || (!isUpper && i < j)) {
            return 0.0
        }
        return this.mat.data[i * this.mat.stride + j]
    }

    public func dims(): (Int64, Int64) {
        return (this.mat.n, this.mat.n)
    }

    // T performs an implicit transpose by returning the receiver inside a Transpose.
    public func t(): Matrix {
        return Transpose(this)
    }

    public func rawTriangular(): Triangular {
        return this.mat
    }

    public func isUpper(): Bool {
        return isUpperUplo(this.mat.uplo)
    }
}

// A RawTriangular can return a blas64.Triangular representation of the receiver.
// Changes to the blas64.Triangular.Data slice will be reflected in the original
// matrix, changes to the N, Stride, Uplo and Diag fields will not.
public interface RawTriangular {
    func rawTriangular(): Triangular
}

public func isUpperUplo(u: Uplo): Bool {
    /*
    case 里只能是常量
     */
    if(u==Upper){
        return true
    }
    else if(u==Lower){
        return false
    }
    else{
        throw Exception(badTriangle)
    }
}

// NewTriDense creates a new Triangular matrix with n rows and columns. If data == nil,
// a new slice is allocated for the backing slice. If len(data) == n*n, data is
// used as the backing slice, and changes to the elements of the returned TriDense
// will be reflected in data. If neither of these is true, NewTriDense will panic.
// NewTriDense will panic if n is zero.
//
// The data must be arranged in row-major order, i.e. the (i*c + j)-th
// element in the data slice is the {i, j}-th element in the matrix.
// Only the values in the triangular portion corresponding to kind are used.
public func newTriDense(n: Int64, kind: TriKind, data: Array<Float64>): TriDense {
    if (n <= 0) {
        if (n == 0) {
            throw Exception(errZeroLength)
        }
        throw Exception("mat: negative dimension")
    }
    if (data.size != 0 && data.size != n * n) {
        throw Exception(errShape)
    }
    if (data.size == 0) {
        let data_ = Array<Float64>(n * n, repeat: 0.0)
        var uplo = Lower
        if (kind == mUpper) {
            uplo = Upper
        }
        return TriDense(
            Triangular(
                n: n,
                stride: n,
                data: data_,
                uplo: uplo,
                diag: NonUnitDiag,
            ),
            n,
        )
    }
    var uplo = Lower
    if (kind == mUpper) {
        uplo = Upper
    }
    return TriDense(
        Triangular(
            n: n,
            stride: n,
            data: data,
            uplo: uplo,
            diag: NonUnitDiag,
        ),
        n,
    )
}
