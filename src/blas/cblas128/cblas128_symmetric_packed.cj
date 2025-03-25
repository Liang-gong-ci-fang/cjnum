package cjnum.blas.cblas128

import cjnum.complex.*
import cjnum.blas.*

public abstract class SymmetricPackedCore {
    public SymmetricPackedCore(
        public let uplo: Uplo,
        public let n: Int64,
        public let data: Array<Complex128>
    ) {}
}

public class SymmetricPacked <: SymmetricPackedCore {
    public SymmetricPacked(
        uplo!: Uplo = 0,
        n!: Int64 = 0,
        data!: Array<Complex128> = []
    ) {
        super(uplo, n, data)
    }
}