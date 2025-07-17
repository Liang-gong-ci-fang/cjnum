package cjnum.mat

import cjnum.blas.blas64.*

/*
因为BasicMatrix只实现了Matrix接口，并不完全实现Dense类型的所有方法，
不能使用仓颉的类型别名实现
所以这里没有使用类型别名，而是重新定义类
*/
public class BasicMatrix<:Matrix{
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
}

public class basicBanded<:Matrix&Banded{
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
    public func tBand(): Banded{ 
        return TransposeBand(this) 
    }
}


