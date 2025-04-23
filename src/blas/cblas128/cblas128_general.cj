package cjnum.blas.cblas128

import cjnum.complex.*

public abstract class GeneralCore {
    public GeneralCore(
        public let rows: Int64,
        public let cols: Int64,
        public let stride: Int64,
        public let data: Array<Complex128>
    ) {}
}

public class General <: GeneralCore {
    public General(
        rows!: Int64 = 0,
        cols!: Int64 = 0,
        stride!: Int64 = 0,
        data!: Array<Complex128> = []
    ) {
        super(rows, cols, stride, data)
    }

    public func from(a: GeneralCols): Unit {
        if (this.rows != a.rows || this.cols != a.cols) {
            throw Exception("cblas128: mismatched dimension")
        }

        if (this.data.size < (this.rows - 1) * this.stride + this.cols) {
            throw Exception("cblas128: short data slice")
        }

        for (j in 0..a.cols) {
            for (i in 0..a.rows) {
                this.data[i * this.stride + j] = a.data[j * a.stride + i]
            }
        }
    }
}

