package cjnum.blas.cblas128

import cjnum.blas.*
import cjnum.complex.*

public class TriangularCols <: TriangularCore {
    public TriangularCols(
        n!: Int64 = 0,
        stride!: Int64 = 0,
        data!: Array<Complex128> = [],
        uplo!: Uplo = 0,
        diag!: Diag = 0
    ) {
        super(n, stride, data, uplo, diag)
    }

    public func from(a: Triangular): Unit {
        if (this.n != a.n) {
            throw Exception("cblas128: mismatched dimension")
        }

        if (this.uplo != a.uplo) {
            throw Exception("cblas128: mismatched BLAS uplo")
        }

        if (this.diag != a.diag) {
            throw Exception("cblas128: mismatched BLAS diag")
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
        } else if (a.uplo == All) {
            for (i in 0..a.n) {
                for (j in 0..a.n) {
                    this.data[i + j * this.stride] = a.data[i * a.stride + j]
                }
            }
        } else {
            throw Exception("cblas128: bad BLAS uplo")
        }
    }
}
