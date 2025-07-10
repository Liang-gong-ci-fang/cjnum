package cjnum.mat

import cjnum.blas.*
import cjnum.blas.blas64.*

// A RawTriBander can return a blas64.TriangularBand representation of the receiver.
// Changes to the blas64.TriangularBand.Data slice will be reflected in the original
// matrix, changes to the N, K, Stride, Uplo and Diag fields will not.
public interface RawTriBander {
    func rawTriBand(): TriangularBand
}

// TriBandDense represents a triangular band matrix in dense storage format.
public class TriBandDense <: Matrix & RawTriBander {
    let mat: TriangularBand
    public init(mat: TriangularBand) {
        this.mat = mat
    }

    @When[feature != "bounds"]
    public func at(i: Int64, j: Int64): Float64 {
        if (i < 0 || UInt64(i) >= UInt64(this.mat.N)) {
            throw Exception(ErrRowAccess)
        }
        if (j < 0 || UInt64(j) >= UInt64(this.mat.N)) {
            throw Exception(ErrColAccess)
        }
        return this._at(i, j)
    }

    @When[feature != "bounds"]
    public func _at(i: Int64, j: Int64): Float64 {
        let isUpper = this.isUpper()
        if ((isUpper && i > j) || (!isUpper && i < j)) {
            return 0.0
        }
        var kl = this.mat.k
        var ku = 0
        if (isUpper) {
            ku = this.mat.k
            kl = 0
        }
        let pj = j + kl - i
        if (pj < 0 || kl + ku + 1 <= pj) {
            return 0.0
        }
        return this.mat.data[i * this.mat.stride + pj]
    }

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
        let isUpper = this.isUpper()
        if ((isUpper && i > j) || (!isUpper && i < j)) {
            return 0.0
        }
        var kl = this.mat.k
        var ku = 0
        if (isUpper) {
            kl = 0
            ku = this.mat.k
            let pj = j + kl - i
            if (pj < 0 || kl + ku + 1 <= pj) {
                return 0.0
            }
            return this.mat.data[i * this.mat.stride + pj]
        }
        let pj = j + kl - i
        if (pj < 0 || kl + ku + 1 <= pj) {
            return 0.0
        }
        return this.mat.data[i * this.mat.stride + pj]
    }

    // Dims returns the number of rows and columns in the matrix.
    public func dims(): (Int64, Int64) {
        return (this.mat.n, this.mat.n)
    }

    // T performs an implicit transpose by returning the receiver inside a Transpose.
    public func t(): Matrix {
        return Transpose(this)
    }

    public func isUpper(): Bool {
        return isUpperUplo(this.mat.uplo)
    }

    // RawTriBand returns the underlying blas64.TriangularBand used by the receiver.
    // Changes to the blas64.TriangularBand.Data slice will be reflected in the original
    // matrix, changes to the N, K, Stride, Uplo and Diag fields will not.
    public func rawTriBand(): TriangularBand {
        return this.mat
    }
}

// NewTriBandDense creates a new triangular banded matrix with n rows and columns,
// k bands in the direction of the specified kind. If data == nil,
// a new slice is allocated for the backing slice. If len(data) == n*(k+1),
// data is used as the backing slice, and changes to the elements of the returned
// TriBandDense will be reflected in data. If neither of these is true, NewTriBandDense
// will panic. k must be at least zero and less than n, otherwise NewTriBandDense will panic.
//
// The data must be arranged in row-major order constructed by removing the zeros
// from the rows outside the band and aligning the diagonals. For example, if
// the upper-triangular banded matrix
//
//	1  2  3  0  0  0
//	0  4  5  6  0  0
//	0  0  7  8  9  0
//	0  0  0 10 11 12
//	0  0  0 0  13 14
//	0  0  0 0  0  15
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
// which is passed to NewTriBandDense as []float64{1, 2, ..., 15, *, *, *}
// with k=2 and kind = mat.Upper.
// The lower triangular banded matrix
//
//	1  0  0  0  0  0
//	2  3  0  0  0  0
//	4  5  6  0  0  0
//	0  7  8  9  0  0
//	0  0 10 11 12  0
//	0  0  0 13 14 15
//
// becomes (* entries are never accessed)
//   - *  1
//   - 2  3
//     4  5  6
//     7  8  9
//     10 11 12
//     13 14 15
//
// which is passed to NewTriBandDense as []float64{*, *, *, 1, 2, ..., 15}
// with k=2 and kind = mat.Lower.
// Only the values in the band portion of the matrix are used.
public func newTriBandDense(n: Int64, k: Int64, kind: TriKind, data: Array<Float64>): TriBandDense {
    if (n <= 0 || k < 0) {
        if (n == 0) {
            throw Exception(errZeroLength)
        }
        throw Exception(errNegativeDimension)
    }
    if (k + 1 > n) {
        throw Exception(errBandwidth)
    }
    let bc = k + 1
    if (data.size != 0 && data.size != n * bc) {
        throw Exception(errShape)
    }
    if (data.size == 0) {
        let data_ = Array<Float64>(n * bc, item: 0.0)
        var uplo = Lower
        if (kind) {
            uplo = Upper
        }
        return TriBandDense(
            TriangularBand(
                uplo: uplo,
                diag: NonUnitDiag,
                n: n,
                k: k,
                data: data_,
                stride: bc,
            ),
        )
    }
    var uplo = Lower
    if (kind) {
        uplo = Upper
    }
    return TriBandDense(
        TriangularBand(
            uplo: uplo,
            diag: NonUnitDiag,
            n: n,
            k: k,
            data: data,
            stride: bc,
        ),
    )
}
