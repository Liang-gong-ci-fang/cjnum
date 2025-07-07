package cjnum.mat

import cjnum.blas.blas64.*
import cjnum.blas.*
import std.math.*

/*
因为以下注释代码本质上是为了在编译过程中检查是DiagDense类否实现了以下接口，
仓颉编译过程中会对未实现的接口方法进行提示，
所以这里没有实现接口判断
*/
// TODO var (
// 	diagDense *DiagDense
// 	_         Matrix          = diagDense
// 	_         allMatrix       = diagDense
// 	_         denseMatrix     = diagDense
// 	_         Diagonal        = diagDense
// 	_         MutableDiagonal = diagDense
// 	_         Triangular      = diagDense
// 	_         TriBanded       = diagDense
// 	_         Symmetric       = diagDense
// 	_         SymBanded       = diagDense
// 	_         Banded          = diagDense
// 	_         RawBander       = diagDense
// 	_         RawSymBander    = diagDense

// 	diag Diagonal
// 	_    Matrix     = diag
// 	_    Diagonal   = diag
// 	_    Triangular = diag
// 	_    TriBanded  = diag
// 	_    Symmetric  = diag
// 	_    SymBanded  = diag
// 	_    Banded     = diag
// )

// Diagonal represents a diagonal matrix, that is a square matrix that only
// has non-zero terms on the diagonal.
public interface Diagonal {
    func at(i: Int64, j: Int64): Float64
    func dims(): (Int64, Int64)
    func t(): Matrix
    // Diag returns the number of rows/columns in the matrix.
    func diag(): Int

    // The following interfaces are included in the Diagonal
    // interface to allow the use of Diagonal types in
    // functions operating on these types.
    // Banded
    // SymBanded
    // Symmetric
    // Triangular
    // TriBanded
}

// DiagDense represents a diagonal matrix in dense storage format.
public class DiagDense <: Matrix & Diagonal {
    var mat: Vector
    public init(mat: Vector) {
        this.mat = mat
    }

    // Diag returns the dimension of the receiver.
    public func diag(): Int64 {
        return this.mat.n
    }

    // Dims returns the dimensions of the matrix.
    public func dims(): (Int64, Int64) {
        return (this.mat.n, this.mat.n)
    }

    // T returns the transpose of the matrix.
    public func t(): Matrix {
        return this
    }

    // At returns the element at row i, column j.
    @When[feature != "bounds"]
    public func at(i: Int64, j: Int64): Float64 {
        if (i < 0 || UInt64(i) >= UInt64(this.mat.n)) {
            panic(errRowAccess)
        }
        if (j < 0 || UInt64(j) >= UInt64(this.mat.n)) {
            panic(errColAccess)
        }
        return this.at(i, j)
    }

    @When[feature != "bounds"]
    public func _at(i: Int64, j: Int64): Float64 {
        if (i != j) {
            return 0.0
        }
        return this.mat.data[i * this.mat.inc]
    }

    // At returns the element at row i, column j.
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
        if (i != j) {
            return 0.0
        }
        return this.mat.data[i * this.mat.inc]
    }

    // Zero sets all of the matrix elements to zero.
    public func zero() {
        for (i in 0..this.mat.n) {
            this.mat.data[this.mat.inc * i] = 0.0
        }
    }

    // SetDiag sets the element at row i, column i to the value v.
    // It panics if the location is outside the appropriate region of the matrix.
    @When[feature!="bounds"]
    public func setDiag(i: Int64, v: Float64):Unit {
        if (i<0||UInt64(i) >= UInt64(this.mat.n)) {
            throw Exception(errRowAccess)
        }
        d._setDiag(i, v)
    }

    @When[feature!="bounds"]
    public func _setDiag(i: Int64, v: Float64):Unit {
        this.mat.data[i * this.mat.inc] = v
    }

    // SetDiag sets the element at row i, column i to the value v.
    // It panics if the location is outside the appropriate region of the matrix.
    @When[feature=="bounds"]
    public func setDiag(i: Int64, v: Float64):Unit {
        this._setDiag(i, v)
    }

    @When[feature=="bounds"]
    public func _setDiag(i: Int64, v: Float64):Unit {
        if (i<0||UInt64(i) >= UInt64(this.mat.n)) {
            throw Exception(errRowAccess)
        }
        this.mat.data[i * this.mat.inc] = v
    }

    // DiagFrom copies the diagonal of m into the receiver. The receiver must
    // be min(r, c) long or empty, otherwise DiagFrom will panic.
    public func diagFrom(m: Matrix): Unit {
        let (m1, m2) = m.dims()
        let n = min(m1, m2)
        this.reuseAsNonZeroed(n)

        var vec = Vector()
        if (m is DiagDense) {
            let r = m as DiagDense
            vec = r.getOrThrow().mat
        } else if (m is RawBander) {
            let r = m as RawBander
            let mat = r.getOrThrow().rawBand()
            vec = Vector(
                n: n,
                inc: mat.stride,
                data: mat.data[mat.kl..(n - 1) * mat.stride + mat.kl + 1],
            )
        } else if (m is RawMatrixer) {
            let r = m as RawMatrixer
            let mat = r.getOrThrow().rawMatrix()
            vec = Vector(
                n: n,
                inc: mat.stride + 1,
                data: mat.data[0..(n - 1) * mat.stride + n],
            )
        } else if (m is RawSymBander) {
            let r = m as RawSymBander
            let mat = r.getOrThrow().rawSymBand()
            vec = Vector(
                n: n,
                inc: mat.stride + 1,
                data: mat.data[0..(n - 1) * mat.stride + n],
            )
        } else if (m is RawSymmetricer) {
            let r = m as RawSymmetricer
            let mat = r.getOrThrow().rawSymmetric()
            vec = Vector(
                n: n,
                inc: mat.stride + 1,
                data: mat.data[0..(n - 1) * mat.stride + n],
            )
        } else if (m is RawTriBander) {
            let r = m as RawTriBander
            let mat = r.getOrThrow().rawTriBand()
            var data = mat.data
            if (mat.uplo == Lower) {
                data = data[mat.k..]
            }
            vec = Vector(
                n: n,
                inc: mat.stride,
                data: data[..(n - 1) * mat.stride + 1],
            )
        } else if (m is RawTriangular) {
            let r = m as RawTriangular
            let mat = r.getOrThrow().rawTriangular()
            if (mat.diag == UnitDiag) {
                var i = 0
                while (i < n) {
                    this.mat.data[i] = 1.0
                    i += this.mat.inc
                }
                return
            }
            vec = Vector(
                n: n,
                inc: mat.stride + 1,
                data: mat.data[0..(n - 1) * mat.stride + n],
            )
        } else if (m is RawVectorer) {
            let r = m as RawVectorer
            this.mat.data[0] = r.getOrThrow().rawVector().data[0]
            return
        } else {
            for (i in 0..n) {
                this.setDiag(i, m.at(i, i))
            }
            return
        }
        copy(vec, this.mat)
    }

    // reuseAsNonZeroed resizes an empty diagonal to a r×r diagonal,
    // or checks that a non-empty matrix is r×r.
    public func reuseAsNonZeroed(r: Int64): Unit {
        if (r == 0) {
            throw Exception(errZeroLength)
        }
        if (this.isEmpty()) {
            this.mat = Vector(
                n: 0,
                inc: 1,
                data: use(this.mat.data, r),
            )
            this.mat.n = r
            return
        }
        if (r != this.mat.n) {
            throw Exception(errShape)
        }
    }

    // IsEmpty returns whether the receiver is empty. Empty matrices can be the
    // receiver for size-restricted operations. The receiver can be emptied using
    // Reset.
    public func isEmpty(): Bool {
        // It must be the case that d.Dims() returns
        // zeros in this case. See comment in Reset().
        return this.mat.inc == 0
    }
}

// NewDiagDense creates a new Diagonal matrix with n rows and n columns.
// The length of data must be n or data must be nil, otherwise NewDiagDense
// will panic. NewDiagDense will panic if n is zero.
public func newDiagDense(n: Int64, data_: Array<Float64>): DiagDense {
    if (n <= 0) {
        if (n == 0) {
            throw Exception(errZeroLength)
        }
        throw Exception("mat: negative dimension")
    }
    if (data_.size == 0) {
        var data = Array<Float64>(n, item: 0.0)
        return DiagDense(Vector(n: n, data: data, inc: 1),)
    }
    if (data_.size != n) {
        throw Exception(errShape)
    }
    return DiagDense(Vector(n: n, data: data_, inc: 1),)
}

/*
原位置在mat/diagonal_test.cj，因为在仓颉中，文件名以 _test.cj 结尾的文件在
正常构建时，二进制文件或库的构建不包含这些文件,
所以这里把DiagDenseViewer类定义的位置放到mat/diagonal.cj中
*/
// diagDenseViewer takes the view of the Diagonal with the underlying Diagonal
// as the DiagDense type.
public interface DiagDenseViewer <: Matrix {
    func at(i: Int64, j: Int64): Float64
    func dims(): (Int64, Int64)
    func t(): Matrix
    func diagView(): Diagonal
}
