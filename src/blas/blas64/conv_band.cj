package cjnum.blas.blas64

public class BandCols <: BandCore {
    public BandCols(
        rows!: Int64 = 0,
        cols!: Int64 = 0,
        kl!: Int64 = 0,
        ku!: Int64 = 0,
        stride!: Int64 = 0,
        data!: Array<Float64> = []
    ) {
        super(rows, cols, kl, ku, stride, data)
    }

    public func from(a: Band): Unit {
        if (this.rows != a.rows || this.cols != a.cols) {
            throw Exception("blas64: mismatched dimension")
        }

        if (this.kl != a.kl || this.ku != a.ku) {
            throw Exception("blas64: mismatched bandwidth")
        }

        if (a.stride < a.kl + a.ku + 1) {
            throw Exception("blas64: short stride for source")
        }

        if (this.stride < this.kl + this.ku + 1) {
            throw Exception("blas64: short stride for destination")
        }

        for (i in 0..a.rows) {
            for (j in max(0, i - a.kl)..min(i + a.ku + 1, a.cols)) {
                this.data[i + this.ku - j + j * this.stride] = a.data[j + a.kl - i + i * a.stride]   
            }
        }
    }
}
