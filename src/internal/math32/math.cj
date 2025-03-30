package cjnum.internal.math32

import std.math.*
import std.math.numeric.*

const mask = 0x7f8 >> 3
const shift = 32 - 8 - 1
const bias = 127

public func nAbs(x: Float32): Float32 {
    match {
        case x < 0.0 => return -x
        case x == 0.0 => return 0.0
        case _ => return x
    }
}

public func nCopysign(x: Float32, y: Float32): Float32 {
    const sign: UInt32 = 1 << 31
    let xBits = x.toBits()
    let yBits = y.toBits()
    let resultBits = (xBits & (sign ^ 0xFFFFFFFF)) | (yBits & sign)

    return Float32.fromBits(resultBits)
}

public func nInf(sign: Int): Float32 {
    var v: Float32
    if (sign >= 0) {
        v = Float32.Inf
    } else {
        v = -Float32.Inf
    }
    return v
}

public func nIsInf(f: Float32, sign: Int): Bool {
    return sign >= 0 && f == Float32.Inf || sign <= 0 && f == -Float32.Inf
}

public func nIsNaN(f: Float32): Bool {
    return f != f
}

public func nHypot(p: Float32, q: Float32): Float32 {
    let x: Float32 = Float32.NaN

    if (nIsInf(p, 0) || nIsInf(q, 0)) {
        return Float32.Inf
    }
    if (nIsNaN(p) || nIsNaN(q)) {
        return Float32.NaN
    }

    var m = nAbs(p)
    var n = nAbs(q)
    if (m < n) {
        let temp = m;
        m = n
        n = temp
    }
    if (m == 0.0) {
        return 0.0
    }
    let result = n / m
    return m * sqrt(1.0 + result * result)
}

public func nSignbit(x: Float32): Bool {
    if (x < 0.0) {
        return true
    }
    if (x == 0.0 && (1.0 / x) < 0.0) {
        return true
    }
    return false
}

public func nMax(x: Float32, y: Float32): Float32 {
    if (nIsInf(x, 1) || nIsInf(y, 1)) {
        return Float32.Inf
    }

    if (nIsNaN(x) || nIsNaN(y)) {
        return Float32.NaN
    }

    if (x == 0.0 && x == y) {
        if (nSignbit(x)) {
            return y
        }
        return x
    }

    if (x > y) {
        return x
    }
    return y
}

public func nMin(x: Float32, y: Float32): Float32 {
    if (nIsInf(x, -1) || nIsInf(y, -1)) {
        return Float32.Inf
    }

    if (nIsNaN(x) || nIsNaN(y)) {
        return Float32.NaN
    }

    if (x == 0.0 && x == y) {
        if (nSignbit(x)) {
            return x
        }
        return y
    }

    if (x < y) {
        return x
    }
    return y
}

public func nNaN(): Float32 {
    return Float32.NaN
}
