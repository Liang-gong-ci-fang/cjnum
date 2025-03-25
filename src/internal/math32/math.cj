package cjnum.internal.math32

import std.math.*
import std.math.numeric.*

const unan = 0x7fc00000
const uinf: UInt32 = 0x7f800000
const uneginf: UInt32 = 0xff800000
const mask = 0x7f8 >> 3
const shift = 32 - 8 - 1
const bias = 127

func Abs(x: Float32): Float32 {
    match {
        case x < 0.0 => return -x
        case x == 0.0 => return 0.0
        case _ => return x
    }
}

func Copysign(x: Float32, y: Float32): Float32 {
    const sign: UInt32 = 1 << 31
    let xBits = x.toBits()
    let yBits = y.toBits()

    let resultBits = (xBits & (sign ^ 0xFFFFFFFF)) | (yBits & sign)

    return Float32.fromBits(resultBits)
}

func Inf(sign: Int): Float32 {
    var v: UInt32
    if (sign >= 0) {
        v = uinf
    } else {
        v = uneginf
    }
    return Float32.fromBits(v)
}

func IsInf(f: Float32, sign: Int): Bool {
    return sign >= 0 && f == Float32(uinf) || sign <= 0 && f == Float32(uneginf)
}

func IsNaN(f: Float32): Bool {
    return f != f
}

func Hypot(p: Float32, q: Float32): Float32 {
    let x: Float32 = Float32.NaN
    
    if (IsInf(p, 0) || IsInf(q, 0)) {
        return Float32(uinf)
    }
    if (IsNaN(p) || IsNaN(q)) {
        return Float32(unan)
    }
    var m = Abs(p)
    var n = Abs(q)
    if (m < n) {
        let temp = m;
        m = n;
        n = temp;
    }
    if (m == 0.0) {
        return 0.0
    }
    let result = n / m
    return m * sqrt(1.0 + result * result)
}

func Signbit(x: Float32): Bool {
    if (x < 0.0) {
        return true
    }
    if (x == 0.0 && (1.0 / x) < 0.0) {
        return true
    }
    return false
}

func Max(x: Float32, y: Float32): Float32 {
    if (IsInf(x, 1) || IsInf(y, 1)) {
        return Float32(uinf)
    }

    if (IsNaN(x) || IsNaN(y)) {
        return Float32(unan)
    }

    if (x == 0.0 && x == y) {
        if (Signbit(x)) {
            return y
        }
        return x
    }

    if (x > y) {
        return x
    }
    return y
}

func Min(x: Float32, y: Float32): Float32 {
    if (IsInf(x, -1) || IsInf(y, -1)) {
        return Float32(uinf)
    }

    if (IsNaN(x) || IsNaN(y)) {
        return Float32(unan)
    }

    if (x == 0.0 && x == y) {
        if (Signbit(x)) {
            return x
        }
        return y
    }

    if (x < y) {
        return x
    }
    return y
}

func NaN(): Float32 {
    return Float32(unan)
}
