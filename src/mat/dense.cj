package cjnum.mat

import cjnum.blas.blas64.*

/*
因为以下注释代码本质上是为了在编译过程中检查是Dense类否实现了以下接口，
仓颉编译过程中会对未实现的接口方法进行提示，
所以这里没有实现接口判断
 */
// TODO var (
// 	dense *Dense

// 	_ Matrix      = dense
// 	_ allMatrix   = dense
// 	_ denseMatrix = dense
// 	_ Mutable     = dense

// 	_ ClonerFrom   = dense
// 	_ RowViewer    = dense
// 	_ ColViewer    = dense
// 	_ RawRowViewer = dense
// 	_ Grower       = dense

// 	_ RawMatrixSetter = dense
// 	_ RawMatrixer     = dense

// 	_ Reseter = dense
// )

// Dense is a dense matrix representation.
public class Dense <: Matrix & RawMatrixer {
    let mat: General
    let capRows: Int64
    let capCols: Int64
    public init(mat: General, capRows: Int64, capCols: Int64) {
        this.mat = mat
        this.capRows = capRows
        this.capCols = capCols
    }

    // Dims returns the number of rows and columns in the matrix.
    public func dims(): (Int64, Int64) {
        return (this.mat.rows, this.mat.cols)
    }

    /*
    func (m *Dense) At 存在两种At方法
     */
    // At returns the element at row i, column j.
    @When[feature != "bounds"]
    public func at(i: Int64, j: Int64): Float64 {
        if (i < 0 || UInt64(i) >= UInt64(m.mat.Rows)) {
            throw Exception(errRowAccess)
        }
        if (j < 0 || UInt64(j) >= UInt64(m.mat.Cols)) {
            throw Exception(errColAccess)
        }
        return m._at(i, j)
    }

    @When[feature != "bounds"]
    public func _at(i: Int64, j: Int64): Float64 {
        return m.mat.Data[i * m.mat.Stride + j]
    }

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
        return this.mat.data[i * this.mat.stride + j]
    }

    // T performs an implicit transpose by returning the receiver inside a Transpose.
    public func t(): Matrix {
        return Transpose(this)
    }

    // RawMatrix returns the underlying blas64.General used by the receiver.
    // Changes to elements in the receiver following the call will be reflected
    // in returned blas64.General.
    public func rawMatrix(): General {
        return this.mat
    }
}

// NewDense creates a new Dense matrix with r rows and c columns. If data == nil,
// a new slice is allocated for the backing slice. If len(data) == r*c, data is
// used as the backing slice, and changes to the elements of the returned Dense
// will be reflected in data. If neither of these is true, NewDense will panic.
// NewDense will panic if either r or c is zero.
//
// The data must be arranged in row-major order, i.e. the (i*c + j)-th
// element in the data slice is the {i, j}-th element in the matrix.
public func newDense(r: Int64, c: Int64, data_: Array<Float64>): Dense {
    if (r <= 0 || c <= 0) {
        if (r == 0 || c == 0) {
            throw Exception(errZeroLength)
        }
        throw Exception(errNegativeDimension)
    }
    var data = data_
    if (!data.isEmpty() && r * c != data.size) {
        throw Exception(errShape)
    }
    if (data.isEmpty()) {
        data = Array<Float64>(r * c, repeat: 0.0)
    }
    return Dense(
        General(
            rows: r,
            cols: c,
            stride: c,
            data: data,
        ),
        r,
        c,
    )
}
