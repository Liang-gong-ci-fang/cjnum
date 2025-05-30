package cjnum.internal.asm.f64

public func dotUnitary(x: Array<Float64>, y: Array<Float64>): Float64 {
    var sum: Float64 = 0.0
    for (i in 0..x.size) {
        sum += y[i] * x[i]
    }
    return sum
}

public func dotInc(x: Array<Float64>, y: Array<Float64>, n: Int64, incX: Int64, incY: Int64, ix: Int64, iy: Int64): Float64 {
    var sum: Float64 = 0.0
    var ix_ = ix
    var iy_ = iy
    for (_ in 0..n) {
        sum += y[iy_] * x[ix_]
        ix_ += incX
        iy_ += incY
    }
    return sum
}
