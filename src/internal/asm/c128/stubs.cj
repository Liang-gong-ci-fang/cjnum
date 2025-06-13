package cjnum.internal.asm.c128

import cjnum.complex.*
import std.math.*
import std.math.numeric.*

public func add(dst: Array<Complex128>, s: Array<Complex128>) {
    for (i in 0..s.size) {
        dst[i] += s[i]
    }
}

public func addConst(alpha: Complex128, x: Array<Complex128>) {
    for (i in 0..x.size) {
        x[i] += alpha
    }
}

public func cumSum(dst: Array<Complex128>, s: Array<Complex128>): Array<Complex128> {
    if (s.size == 0) {
        return dst
    }
    dst[0] = s[0]
    for (i in 1..s.size) {
        dst[i] = dst[i - 1] + s[i]
    }
    return dst
}

public func cumProd(dst: Array<Complex128>, s: Array<Complex128>): Array<Complex128> {
    if (s.size == 0) {
        return dst
    }
    dst[0] = s[0]
    for (i in 1..s.size) {
        dst[i] = dst[i - 1] * s[i]
    }
    return dst
}

public func div(dst: Array<Complex128>, s: Array<Complex128>) {
    for (i in 0..s.size) {
        dst[i] /= s[i]
    }
}

public func divTo(dst: Array<Complex128>, s: Array<Complex128>, t: Array<Complex128>): Array<Complex128> {
    for (i in 0..s.size) {
        dst[i] = s[i] / t[i]
    }
    return dst
}

public func dotUnitary(x: Array<Complex128>, y: Array<Complex128>): Complex128 {
    var sum: Complex128 = Complex128(0.0)
    for (i in 0..x.size) {
        sum += x[i].conj() * y[i]
    }
    return sum
}

public func l2DistanceUnitary(x: Array<Complex128>, y: Array<Complex128>): Float64 {
    var scale: Float64 = 0.0
    var sumSquares: Float64 = 1.0
    for (i in 0..x.size) {
        var v = x[i] - y[i]
        if (v == Complex128(0.0)) {
            continue
        }
        let absxi = c128Abs(v)
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

public func l2NormUnitary(x: Array<Complex128>): Float64 {
    var scale: Float64 = 0.0
    var sumSquares: Float64 = 1.0
    var s: Float64 = 0.0
    for (v in x) {
        if (v == Complex128(0.0)) {
            continue
        }
        let absxi = c128Abs(v)
        if (absxi.isNaN()) {
            return Float64.NaN
        }
        if (scale < absxi) {
            s = scale / absxi
            sumSquares = 1.0 + sumSquares * s * s
            scale = absxi
        } else {
            s = absxi / scale
            sumSquares += s * s
        }
    }
    if (scale.isInf()) {
        return Float64.Inf
    }
    return scale * sqrt(sumSquares)
}

public func sum(x: Array<Complex128>): Complex128 {
    var sum: Complex128 = Complex128(0.0)
    for (v in x) {
        sum += v
    }
    return sum
}
