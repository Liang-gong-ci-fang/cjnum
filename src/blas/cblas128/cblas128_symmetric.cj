package cjnum.blas.cblas128

import cjnum.complex.*
import cjnum.blas.*

public abstract class SymmetricCore {
    public SymmetricCore(
        public let n: Int64,
        public let stride: Int64,
        public let data: Array<Complex128>,
        public var uplo: Uplo
    ) {}
}

public class Symmetric <: SymmetricCore {
    public Symmetric(
        n!: Int64 = 0,
        stride!: Int64 = 0,
        data!: Array<Complex128> = [],
        uplo!: Uplo = 0
    ) {
        super(n, stride, data, uplo)
    }

    public func from(a: SymmetricCols): Unit {
        if (this.n != a.n) {
            throw Exception("cblas128: mismatched dimension")
        }
        if (this.uplo != a.uplo) {
            throw Exception("cblas128: mismatched BLAS uplo")
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
        } else {
            throw Exception("cblas128: bad BLAS uplo")
        }
    }
}