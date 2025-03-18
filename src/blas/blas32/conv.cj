package cjnum.blas.blas32

public class GeneralCols <: GeneralCore {
    public GeneralCols(
        rows: Int64,
        cols: Int64,
        stride: Int64,
        data: Array<Float32>
    ) {
        super(rows, cols, stride, data)
    }

    public func from(a: General): Unit {
        if (this.rows != a.rows || this.cols != a.cols) {
            throw Exception("blas32: mismatched dimension")
        }

        if (this.data.size < (this.cols - 1) * this.stride + this.rows) {
            throw Exception("blas32: short data slice")
        }

        for (i in 0..a.rows) {
            for (j in 0..a.cols) {
                this.data[i + j * this.stride] = a.data[i * a.stride + j]
            }
        }
    }
}
