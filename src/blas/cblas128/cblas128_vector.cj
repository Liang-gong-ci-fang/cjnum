package cjnum.blas.cblas128

import cjnum.complex.*

public abstract class VectorCore {
    public VectorCore(
        public let n: Int64,
        public let data: Array<Complex128>,
        public let inc: Int64
    ) {}
}

public class Vector <: VectorCore {
    public Vector(
        n!: Int64 = 0,
        data!: Array<Complex128> = [],
        inc!: Int64 = 0
    ) {
        super(n, data, inc)
    }
}
