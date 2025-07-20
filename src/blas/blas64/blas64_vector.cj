package cjnum.blas.blas64

import cjnum.blas.*

public abstract class VectorCore {
    public VectorCore(
        public var n: Int64,
        public var data: Array<Float64>,
        public let inc: Int64
    ) {}
}

public class Vector <: VectorCore {
    public Vector(
        n!: Int64 = 0,
        data!: Array<Float64> = [],
        inc!: Int64 = 0
    ) {
        super(n, data, inc)
    }
}