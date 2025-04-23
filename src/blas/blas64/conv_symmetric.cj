package cjnum.blas.blas64

import cjnum.blas.*

public class SymmetricCols <: SymmetricCore {
    public SymmetricCols(
        n!: Int64 = 0,
        stride!: Int64 = 0,
        data!: Array<Float64> = [],
        uplo!: Uplo = 0
    ) {
        super(n, stride, data, uplo)
    }

    public func from(a: Symmetric): Unit {
        if (this.n != a.n) {
            throw Exception("blas64: mismatched dimension")
        }
        if (this.uplo != a.uplo) {
            throw Exception("blas64: mismatched BLAS uplo")
        }
        if (a.uplo == Upper) {
            for (i in 0..a.n) {
                for (j in i..a.n) {
                    this.data[i + j * this.stride] = a.data[i * a.stride + j]
                }
            }
        } else if (a.uplo == Lower) {
            for (i in 0..a.n) {
                for (j in 0..=i) {
                    this.data[i + j * this.stride] = a.data[i * a.stride + j]
                }
            }
        } else {
            throw Exception("blas64: bad BLAS uplo")
        }
    }
}