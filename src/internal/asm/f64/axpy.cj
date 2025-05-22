package cjnum.internal.asm.f64

public func axpyUnitary(alpha: Float64, x: Array<Float64>, y: Array<Float64>): Unit {
    for (i in 0..x.size) {
        y[i] += alpha * x[i]
    }
}

public func axpyUnitaryTo(dst: Array<Float64>, alpha: Float64, x: Array<Float64>, y: Array<Float64>): Unit {
    for (i in 0..x.size) {
        dst[i] = alpha * x[i] + y[i]
    }
}

public func axpyInc(alpha: Float64, x: Array<Float64>, y: Array<Float64>, n: Int64, incX: Int64, incY: Int64, ix: Int64,
    iy: Int64): Unit {
    var ix_ = ix
    var iy_ = iy
    for (_ in 0..n) {
        y[iy_] += alpha * x[ix_]
        ix_ += incX
        iy_ += incY
    }
}

public func axpyIncTo(dst: Array<Float64>, incDst: Int64, idst: Int64, alpha: Float64, x: Array<Float64>,
    y: Array<Float64>, n: Int64, incX: Int64, incY: Int64, ix: Int64, iy: Int64): Unit {
    var ix_ = ix
    var iy_ = iy
    var idst_ = idst
    for (_ in 0..n) {
        dst[idst_] = alpha * x[ix_] + y[iy_]
        ix_ += incX
        iy_ += incY
        idst_ += incDst
    }
}
