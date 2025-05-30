package cjnum.internal.asm.f64

import std.math.*
import std.math.numeric.*

public func l2NormUnitary(x: Array<Float64>): Float64 {
    var scale: Float64 = 0.0
    var sumSquares: Float64 = 1.0
    for (v in x) {
        if (v == 0.0) {
            continue
        }
        let absxi = abs(v)
        if (absxi.isNaN()) {
            return Float64.NaN
        }
        if (scale < absxi) {
            var s = scale / absxi
            sumSquares = 1.0 + sumSquares * s * s
            scale = absxi
        } else {
            var s = absxi / scale
            sumSquares += s * s
        }
    }
    if (scale == Float64.Inf) {
        return Float64.Inf
    }
    return scale * sqrt(sumSquares)
}

public func l2NormInc(x: Array<Float64>, n: Int64, incX: Int64): Float64 {
    var scale: Float64 = 0.0
    var sumSquares: Float64 = 1.0
    var ix = 0
    while (ix < n * incX) {
        var val = x[ix]
        if (val == 0.0) {
            ix += incX
            continue
        }
        let absxi = abs(val)
        if (absxi.isNaN()) {
            return Float64.NaN
        }
        if (scale < absxi) {
            var s = scale / absxi
            sumSquares = 1.0 + sumSquares * s * s
            scale = absxi
        } else {
            var s = absxi / scale
            sumSquares += s * s
        }
        ix += incX
    }
    if (scale == Float64.Inf) {
        return Float64.Inf
    }
    return scale * sqrt(sumSquares)
}

public func l2DistanceUnitary(x: Array<Float64>, y: Array<Float64>): Float64 {
    var scale: Float64 = 0.0
    var sumSquares: Float64 = 1.0
    for (i in 0..x.size) {
        var v = x[i] - y[i]
        if (v == 0.0) {
            continue
        }
        var absxi = abs(v)
        if (absxi.isNaN()) {
            return Float64.NaN
        }
        if (scale < absxi) {
            var s = scale / absxi
            sumSquares = 1.0 + sumSquares * s * s
            scale = absxi
        } else {
            var s = absxi / scale
            sumSquares += s * s
        }
    }
    if (scale == Float64.Inf) {
        return Float64.Inf
    }
    return scale * sqrt(sumSquares)
}
