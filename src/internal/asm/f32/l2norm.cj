package cjnum.internal.asm.f32

import cjnum.internal.math32.*
import std.math.*

public func l2NormUnitary(x: Array<Float32>): Float32 {
    var scale: Float32 = 0.0
    var sumSquares: Float32 = 1.0
    for (i in 0..x.size) {
        let v = x[i]
        if (v == 0.0) {
            continue
        }
        let absxi = nAbs(v)
        if (nIsNaN(absxi)) {
            return nNaN()
        }
        if (scale < absxi) {
            let s = scale / absxi
            sumSquares = 1.0 + sumSquares * s * s
            scale = absxi
        } else {
            let s = absxi / scale
            sumSquares += s * s
        }
    }
    if (nIsInf(scale, 1)) {
        return nInf(1)
    }
    return scale * sqrt(sumSquares)
}

public func l2NormInc(x: Array<Float32>, n: Int64, incX: Int64): Float32 {
    var scale: Float32 = 0.0
    var sumSquares: Float32 = 1.0
    for (ix in 0..n) {
        let val = x[ix * incX]
        if (val == 0.0) {
            continue
        }
        let absxi = nAbs(val)
        if (nIsNaN(absxi)) {
            return nNaN()
        }
        if (scale < absxi) {
            let s = scale / absxi
            sumSquares = 1.0 + sumSquares * s * s
            scale = absxi
        } else {
            let s = absxi / scale
            sumSquares += s * s
        }
    }
    if (nIsInf(scale, 1)) {
        return nInf(1)
    }
    return scale * sqrt(sumSquares)
}

public func l2DistanceUnitary(x: Array<Float32>, y: Array<Float32>): Float32 {
    var scale: Float32 = 0.0
    var sumSquares: Float32 = 1.0
    for (i in 0..x.size) {
        let v = x[i] - y[i]
        if (v == 0.0) {
            continue
        }
        let absxi = nAbs(v)
        if (nIsNaN(absxi)) {
            return nNaN()
        }
        if (scale < absxi) {
            let s = scale / absxi
            sumSquares = 1.0 + sumSquares * s * s
            scale = absxi
        } else {
            let s = absxi / scale
            sumSquares += s * s
        }
    }
    if (nIsInf(scale, 1)) {
        return nInf(1)
    }
    return scale * sqrt(sumSquares)
}
