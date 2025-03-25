package cjnum.blas.cblas128

import cjnum.complex.*

public class GeneralCols <: GeneralCore {
    public GeneralCols(
        rows!: Int64 = 0,
        cols!: Int64 = 0,
        stride!: Int64 = 0,
        data!: Array<Complex128> = []
    ) {
        super(rows, cols, stride, data)
    }

    public func from(a: General): Unit {
        if (this.rows != a.rows || this.cols != a.cols) {
            throw Exception("cblas128: mismatched dimension")
        }

        if (this.data.size < (this.cols - 1) * this.stride + this.rows) {
            throw Exception("cblas128: short data slice")
        }

        for (i in 0..a.rows) {
            for (j in 0..a.cols) {
                this.data[i + j * this.stride] = a.data[i * a.stride + j]
            }
        }
    }
}
