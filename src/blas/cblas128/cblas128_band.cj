package cjnum.blas.cblas128

import cjnum.complex.*

public abstract class BandCore {
    public BandCore(
        public let rows: Int64,
        public let cols: Int64,
        public var kl: Int64,
        public var ku: Int64,
        public let stride: Int64,
        public let data: Array<Complex128>
    ) {}
}

public class Band <: BandCore {
    public Band(
        rows!: Int64 = 0,
        cols!: Int64 = 0,
        kl!: Int64 = 0,
        ku!: Int64 = 0,
        stride!: Int64 = 0,
        data!: Array<Complex128> = []
    ) {
        super(rows, cols, kl, ku, stride, data)
    }

    public func from(a: BandCols): Unit {
        if (this.rows != a.rows || this.cols != a.cols) {
            throw Exception("cblas128: mismatched dimension")
        }
        if (this.kl != a.kl || this.ku != a.ku) {
            throw Exception("cblas128: mismatched bandwidth")
        }
        if (a.stride < a.kl + a.ku + 1) {
            throw Exception("cblas128: short stride for source")
        }
        if (this.stride < this.kl + this.ku + 1) {
            throw Exception("cblas128: short stride for destination")
        }
        for (j in 0..a.cols) {
            for (i in max(0, j - a.ku)..min(j + a.kl + 1, a.rows)) {
                this.data[j + a.kl - i + i * a.stride] = a.data[i + this.ku - j + j * this.stride]
            }
        }
    }
}
