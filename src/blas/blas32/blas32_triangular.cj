package cjnum.blas.blas32

import cjnum.blas.*

public abstract class TriangularCore {
    public TriangularCore(
        public let n: Int64,
        public let stride: Int64,
        public let data: Array<Float32>,
        public var uplo: Uplo,
        public var diag: Diag
    ) {}
}

public class Triangular <: TriangularCore {
    public Triangular(
        n!: Int64 = 0,
        stride!: Int64 = 0,
        data!: Array<Float32> = [],
        uplo!: Uplo = 0,
        diag!: Diag = 0
    ) {
        super(n, stride, data, uplo, diag)
    }

    public func from(a: TriangularCols): Unit {
        if (this.n != a.n) {
            throw Exception("blas32: mismatched dimension")
        }

        if (this.uplo != a.uplo) {
            throw Exception("blas32: mismatched BLAS uplo")
        }

        if (this.diag != a.diag) {
            throw Exception("blas32: mismatched BLAS diag")
        }

        if (a.uplo == Upper) {
            for (i in 0..a.n) {
                for (j in i..a.n) {
                    this.data[i * this.stride + j] = a.data[i + j * a.stride]
                }
            }
        } else if (a.uplo == Lower) {
            for (i in 0..a.n) {
                for (j in 0..=i) {
                    this.data[i * this.stride + j] = a.data[i + j * a.stride]
                }
            }
        } else if (a.uplo == All) {
            for (i in 0..a.n) {
                for (j in 0..a.n) {
                    this.data[i * this.stride + j] = a.data[i + j * a.stride]
                }
            }
        } else {
            throw Exception("blas32: bad BLAS uplo")
        }
    }
}
