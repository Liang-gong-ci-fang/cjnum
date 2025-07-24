package cjnum.mat

import cjnum.blas.blas64.*
import std.math.*
import std.collection.*

// Matrix is the basic matrix interface type.
public interface Matrix {
    // Dims returns the dimensions of a Matrix.
    func dims(): (Int64, Int64)

    // At returns the value of a matrix element at row i, column j.
    // It will panic if i or j are out of bounds for the matrix.
    func at(i: Int64, j: Int64): Float64

    // T returns the transpose of the Matrix. Whether T returns a copy of the
    // underlying data is implementation dependent.
    // This method may be implemented using the Transpose type, which
    // provides an implicit matrix transpose.
    func t(): Matrix
}

public interface allMatrix <: reseter{
    func isEmpty(): Bool
    func zero(): Unit
}

public interface reseter {
    func reset(): Unit
}

/*
因为以下注释代码本质上是为了在编译过程中检查是Transpose类否实现了以下接口，
仓颉编译过程中会对未实现的接口方法进行提示，
所以这里没有实现接口判断
 */
// TODO var (
// 	_ Matrix       = Transpose{}
// 	_ Untransposer = Transpose{}
// )

// Transpose is a type for performing an implicit matrix transpose. It implements
// the Matrix interface, returning values from the transpose of the matrix within.
public class Transpose <: Matrix & Untransposer {
    let matrix: Matrix
    public init(matrix: Matrix) {
        this.matrix = matrix
    }

    // Dims returns the dimensions of the transposed matrix. The number of rows returned
    // is the number of columns in the Matrix field, and the number of columns is
    // the number of rows in the Matrix field.
    public func dims(): (Int64, Int64) {
        let (c, r) = this.matrix.dims()
        return (r, c)
    }

    // At returns the value of the element at row i and column j of the transposed
    // matrix, that is, row j and column i of the Matrix field.
    public func at(i: Int64, j: Int64): Float64 {
        return this.matrix.at(j, i)
    }

    // T performs an implicit transpose by returning the Matrix field.
    public func t(): Matrix {
        return this.matrix
    }

    // Untranspose returns the Matrix field.
    public func mUntranspose(): Matrix {
        return this.matrix
    }
}

// Untransposer is a type that can undo an implicit transpose.
public interface Untransposer {
    // Note: This interface is needed to unify all of the Transpose types. In
    // the mat methods, we need to test if the Matrix has been implicitly
    // transposed. If this is checked by testing for the specific Transpose type
    // then the behavior will be different if the user uses T() or TTri() for a
    // triangular matrix.

    // Untranspose returns the underlying Matrix stored for the implicit transpose.
    func mUntranspose(): Matrix
}

// untranspose untransposes a matrix if applicable. If a is an Untransposer, then
// untranspose returns the underlying matrix and true. If it is not, then it returns
// the input matrix and false.
private func untranspose(a: Matrix): (Matrix, Bool) {
    let ut = a as Untransposer
    if (!ut.isNone()) {
        let mx = ut.getOrThrow()
        return (mx.mUntranspose(), true)
    }
    println("ut isNone")
    return (a, false)
}

// A RawMatrixer can return a blas64.General representation of the receiver. Changes to the blas64.General.Data
// slice will be reflected in the original matrix, changes to the Rows, Cols and Stride fields will not.
public interface RawMatrixer {
    func rawMatrix(): General
}

// A RawVectorer can return a blas64.Vector representation of the receiver. Changes to the blas64.Vector.Data
// slice will be reflected in the original matrix, changes to the Inc field will not.
public interface RawVectorer {
    func rawVector(): Vector
}

// Col copies the elements in the jth column of the matrix into the slice dst.
// The length of the provided slice must equal the number of rows, unless the
// slice is nil in which case a new slice is first allocated.
public func mCol(dst_: Array<Float64>, j: Int64, a: Matrix): Array<Float64> {
    var dst = dst_
    let (r, c) = a.dims()
    println("r= ${r}")
    if (j < 0 || j >= c) {
        throw Exception(errColAccess)
    }
    if (dst.size == 0) {
        dst = Array<Float64>(r, item: 0.0)
        println("1. dst size= ${dst.size}")
    } else {
        if (dst.size != r) {
            throw Exception(errColLength)
        }
    }
    println("success judge")
    let (aU, aTrans) = untranspose(a)
    let rm = aU as RawMatrixer
    if (rm.isSome()) {
        println("isSome begin")
        let m = rm.getOrThrow().rawMatrix()
        //
        if (aTrans) {
            println("copyTo begin")
            m.data.copyTo(dst, j * m.stride, 0, m.cols)
            return dst
        }
        //
        println("copy begin")
        println("m.stride = ${m.stride} m.data[j..].size = ${m.data[j..].size}")
        println("inc = 1 dst.size = ${dst.size}")
        copy(
            Vector(n: r, inc: m.stride, data: m.data[j..]),
            Vector(n: r, inc: 1, data: dst)
        )
        println("success isSome")
        return dst
    }
    for (i in 0..r) {
        dst[i] = a.at(i, j)
    }
    println("success isNone")
    return dst
}

// Equal returns whether the matrices a and b have the same size
// and are element-wise equal.
public func mEqual(a: Matrix, b: Matrix): Bool {
    let (ar, ac) = a.dims()
    let (br, bc) = b.dims()
    if (ar != br || ac != bc) {
        return false
    }
    let (aU, aTrans) = untranspose(a)
    let (bU, bTrans) = untranspose(b)
    let rma = aU as RawMatrixer
    if (!rma.isNone()) {
        let rmb = bU as RawMatrixer
        if (!rmb.isNone()) {
            let ra = rma.getOrThrow().rawMatrix()
            let rb = rmb.getOrThrow().rawMatrix()
            if (aTrans == bTrans) {
                for (i in 0..ra.rows) {
                    for (j in 0..ra.cols) {
                        if (ra.data[i * ra.stride + j] != rb.data[i * rb.stride + j]) {
                            return false
                        }
                    }
                }
                return true
            }
            for (i in 0..ra.rows) {
                for (j in 0..ra.cols) {
                    if (ra.data[i * ra.stride + j] != rb.data[j * rb.stride + i]) {
                        return false
                    }
                }
            }
            return true
        }
    }
    let rma2 = aU as RawSymmetricer
    if (!rma2.isNone()) {
        let rmb = bU as RawSymmetricer
        if (!rmb.isNone()) {
            let ra = rma2.getOrThrow().rawSymmetric()
            let rb = rmb.getOrThrow().rawSymmetric()
            // Symmetric matrices are always upper and equal to their transpose.
            for (i in 0..ra.n) {
                for (j in i..ra.n) {
                    if (ra.data[i * ra.stride + j] != rb.data[i * rb.stride + j]) {
                        return false
                    }
                }
            }
            return true
        }
    }
    let ra = aU as VecDense
    if (!ra.isNone()) {
        let rb = bU as VecDense
        if (!rb.isNone()) {
            // If the raw vectors are the same length they must either both be
            // transposed or both not transposed (or have length 1).
            for (i in 0..ra.getOrThrow().mat.n) {
                let amat = ra.getOrThrow().mat
                let bmat = rb.getOrThrow().mat
                if (amat.data[i * amat.inc] != bmat.data[i * bmat.inc]) {
                    return false
                }
            }
            return true
        }
    }
    for (i in 0..ar) {
        for (j in 0..ac) {
            if (a.at(i, j) != b.at(i, j)) {
                return false
            }
        }
    }
    return true
}

// zero zeros the given slice's elements.
public func _zero(f_: Array<Float64>): Unit {
    var f = f_
    for (i in 0..f.size) {
        f[i] = 0.0
    }
}

// use returns a float64 slice with l elements, using f if it
// has the necessary capacity, otherwise creating a new slice.
public func use(f_: Array<Float64>, l: Int64): Array<Float64> {
    let f: ArrayList<Float64> = ArrayList<Float64>(f_)
    if (l <= f.capacity()) {
        return f.slice(0..l).toArray()
    }
    return Array<Float64>(l, item: 0.0)
}
