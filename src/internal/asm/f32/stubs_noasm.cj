package cjnum.internal.asm.f32

public func axpyUnitary(alpha: Float32, x: Array<Float32>, y: Array<Float32>) {
    for (i in 0..x.size) {
        y[i] += alpha * x[i]
    }
}

public func axpyUnitaryTo(dst: Array<Float32>, alpha: Float32, x: Array<Float32>, y: Array<Float32>) {
    for (i in 0..x.size) {
        dst[i] = alpha * x[i] + y[i]
    }
}

public func axpyInc(alpha: Float32, x: Array<Float32>, y: Array<Float32>, n: Int64, incX: Int64, incY: Int64, ix: Int64,
    iy: Int64) {
    var ix_ = ix
    var iy_ = iy
    for (i in 0..n) {
        y[iy_] += alpha * x[ix_]
        ix_ += incX
        iy_ += incY
    }
}

public func axpyIncTo(dst: Array<Float32>, incDst: Int64, idst: Int64, alpha: Float32, x: Array<Float32>,
    y: Array<Float32>, n: Int64, incX: Int64, incY: Int64, ix: Int64, iy: Int64) {
    var ix_ = ix
    var iy_ = iy
    var idst_ = idst
    for (i in 0..n) {
        dst[idst_] = alpha * x[ix_] + y[iy_]
        ix_ += incX
        iy_ += incY
        idst_ += incDst
    }
}

public func dotUnitary(x: Array<Float32>, y: Array<Float32>): Float32 {
    var sum: Float32 = 0.0
    for (i in 0..x.size) {
        sum += y[i] * x[i]
    }
    return sum
}

public func dotInc(x: Array<Float32>, y: Array<Float32>, n: Int64, incX: Int64, incY: Int64, ix: Int64, iy: Int64): Float32 {
    var ix_ = ix
    var iy_ = iy
    var sum: Float32 = 0.0
    for (i in 0..n) {
        sum += y[iy_] * x[ix_]
        ix_ += incX
        iy_ += incY
    }
    return sum
}

public func ddotUnitary(x: Array<Float32>, y: Array<Float32>): Float64 {
    var sum: Float64 = 0.0
    for (i in 0..x.size) {
        sum += Float64(y[i]) * Float64(x[i])
    }
    return sum
}

public func ddotInc(x: Array<Float32>, y: Array<Float32>, n: Int64, incX: Int64, incY: Int64, ix: Int64, iy: Int64): Float64 {
    var ix_ = ix
    var iy_ = iy
    var sum: Float64 = 0.0
    for (i in 0..n) {
        sum += Float64(y[iy_]) * Float64(x[ix_])
        ix_ += incX
        iy_ += incY
    }
    return sum
}

public func sum(x: Array<Float32>): Float32 {
    var sum: Float32 = 0.0
    for (i in 0..x.size) {
        sum += x[i]
    }
    return sum
}
