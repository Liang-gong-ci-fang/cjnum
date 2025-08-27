package cjnum.lapack.lapack64

import cjnum.lapack.*
import cjnum.lapack.cjnum.*
import cjnum.blas.blas64.*

var lapack64: LFloat64 = LFloat64Implementation()

public func lFloat64Use(l: LFloat64): Unit {
    lapack64 = l
}

public func langb(norm: MatrixNorm, a: Band): Float64 {
    return lapack64.dlangb(norm, a.rows, a.cols, a.kl, a.ku, a.data, max(1, a.stride))
}