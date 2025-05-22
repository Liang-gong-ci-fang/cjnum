package cjnum.internal.asm.f64

public func scalUnitary(alpha: Float64, x: Array<Float64>) {
    for (i in 0..x.size) {
        x[i] *= alpha
    }
}

public func scalUnitaryTo(dst: Array<Float64>, alpha: Float64, x: Array<Float64>) {
    for (i in 0..x.size) {
        dst[i] = alpha * x[i]
    }
}

public func scalInc(alpha: Float64, x: Array<Float64>, n: Int64, incX: Int64) {
    var ix: Int64 = 0
    for (_ in 0..n) {
        x[ix] *= alpha
        ix += incX
    }
}

public func scalIncTo(dst: Array<Float64>, incDst: Int64, alpha: Float64, x: Array<Float64>, n: Int64, incX: Int64) {
    var idst: Int64 = 0
    var ix: Int64 = 0
    for (_ in 0..n) {
        dst[idst] = alpha * x[ix]
        ix += incX
        idst += incDst
    }
}
