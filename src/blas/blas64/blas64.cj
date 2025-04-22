package cjnum.blas.blas64

import cjnum.*
import cjnum.blas.*
import cjnum.blas.cjnum.*

var blas64: NFloat64 = NFloat64Implementation()

public func nFloat64use(b: NFloat64): Unit {
    blas64 = b
}

public func nFloat64Implementation(): NFloat64 {
    return blas64
}

private const negInc    = "blas64: negative vector increment"
private const badLength = "blas64: vector length mismatch"

// TODO