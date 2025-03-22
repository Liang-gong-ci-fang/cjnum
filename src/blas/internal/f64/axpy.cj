package cjnum.blas.internal.f64

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

public func axpyInc(dst: Array<Float64>, alpha: Float64, x: Array<Float64>, y: Array<Float64>,
                        n: Int64, incX: Int64, incY: Int64, ix: Int64, iy: Int64): Unit {
    var iix = ix
    var iiy = iy
    for (i in 0..n) {
        dst[iiy] += alpha * x[iix]
        iix += incX
        iiy += incY
    }
}

public func axpyIncTo(dst: Array<Float64>, incDst: Int64, idst: Int64, alpha: Float64, x: Array<Float64>, y: Array<Float64>,
                        n: Int64, incX: Int64, incY: Int64, ix: Int64, iy: Int64): Unit {
    var iix = ix
    var iiy = iy
    var iidst = idst
    for (i in 0..n) {
        dst[iidst] = alpha * x[iix] + y[iiy]
        iix += incX
        iiy += incY
        iidst += incDst
    }
}
