package cjnum.blas.blas64

import cjnum.blas.*

public abstract class SymmetricPackedCore {
    public SymmetricPackedCore(
        public let uplo: Uplo,
        public let n: Int64,
        public let data: Array<Float64>
    ) {}
}

public class SymmetricPacked <: SymmetricPackedCore {
    public SymmetricPacked(
        uplo!: Uplo = 0,
        n!: Int64 = 0,
        data!: Array<Float64> = []
    ) {
        super(uplo, n, data)
    }
}