package cjnum.blas.cblas128

import cjnum.complex.*
import cjnum.blas.*

public class SymmetricBandCols <: SymmetricBandCore {
    public SymmetricBandCols(
        n!: Int64 = 0,
        k!: Int64 = 0,
        stride!: Int64 = 0,
        data!: Array<Complex128> = [],
        uplo!: Uplo = 0
    ) {
        super(n, k, stride, data, uplo)
    }

    public func from(a: SymmetricBand): Unit {
        if (this.n != a.n) {
            throw Exception("cblas128: mismatched dimension")
        }
        if (this.k != a.k) {
            throw Exception("cblas128: mismatched bandwidth")
        }
        if (a.stride < a.k + 1) {
            throw Exception("cblas128: short stride for source")
        }
        if (this.stride < this.k + 1) {
            throw Exception("cblas128: short stride for destination")
        }
        if (this.uplo != a.uplo) {
            throw Exception("cblas128: mismatched BLAS uplo")
        }
        let dst = BandCols(rows: this.n, cols: this.n, stride: this.stride, data: this.data)
        let src = Band(rows: a.n, cols: a.n, stride: a.stride, data: a.data)
        if (a.uplo == Upper) {
            dst.ku = this.k
            src.ku = a.k
        } else if (a.uplo == Lower) {
            dst.kl = this.k
            src.kl = a.k
        } else {
            throw Exception("cblas128: bad BLAS uplo")
        }
        dst.from(src)
    }
}