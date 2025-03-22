package cjnum.blas.blas64

import cjnum.blas.*

public abstract class TriangularPackedCore {
    public TriangularPackedCore(
        public let uplo: Uplo,
        public let diag: Diag,
        public let n: Int64,
        public let data: Array<Float64>
    ) { }
}

public class TriangularPacked <: TriangularPackedCore {
    public TriangularPacked(
        uplo!: Uplo = 0,
        diag!: Diag = 0,
        n!: Int64 = 0,
        data!: Array<Float64> = []
    ) { 
        super(uplo, diag, n, data)
    }
}