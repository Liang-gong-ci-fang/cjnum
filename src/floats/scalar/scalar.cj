package cjnum.floats.scalar

import std.math.*
import std.math.numeric.*

public func equalWithinAbs(a: Float64, b: Float64, tol: Float64): Bool {
    return a == b || abs(a - b) <= tol
}

public func equalWithinRel(a: Float64, b: Float64, tol: Float64): Bool {
    if (a == b) {
        return true
    }
    let delta = abs(a - b)
    if (delta <= Float64.MinNormal) {
        return delta <= tol * Float64.MinNormal
    }
    return delta / max(abs(a), abs(b)) <= tol
}

public func equalWithinAbsOrRel(a: Float64, b: Float64, absTol: Float64, relTol: Float64): Bool {
    return equalWithinAbs(a, b, absTol) || equalWithinRel(a, b, relTol)
}

public func equalWithinULP(a: Float64, b: Float64, ulp: UInt64): Bool {
    if (a == b) {
        return true
    }
    if (a.isNaN() || b.isNaN()) {
        return false
    }
    if (a >= 0.0 != b >= 0.0) {
        return abs(a).toBits() + abs(b).toBits() <= ulp
    }
    return ulpDiff(abs(a).toBits(), abs(b).toBits()) <= ulp
}

public func ulpDiff(a: UInt64, b: UInt64): UInt64 {
    if (a > b) {
        return a - b
    }
    return b - a
}

const nanBits: UInt64 = 0x7ff8000000000000
const nanMask: UInt64 = 0xfff8000000000000

public func naNWith(payload: UInt64): Float64 {
    return Float64.fromBits(nanBits | (payload & (!nanMask)))
}

public func naNPayload(f: Float64): (UInt64, Bool) {
    let b = f.toBits()
    if ((b & nanBits) != nanBits) {
        return (0, false)
    }
    return (b & (!nanMask), true)
}

public func parseWithNA(s: String, missing: String): (Float64, Float64) {
    if (s == missing) {
        return (0.0, 0.0)
    }
    let value = Decimal.parse(s).toFloat64()
    return (value, 1.0)
}

public func round(x: Float64, prec: Int64): Float64 {
    if (x == 0.0) {
        return 0.0
    }
    if (prec >= 0 && x == trunc(x)) {
        return x
    }
    let pow = Decimal(10).powWithPrecision(prec, 0)
    let intermed = x * pow.toFloat64()
    if (intermed.isInf()) {
        return x
    }
    if (intermed.isNaN()) {
        return intermed
    }

    let x_ = round(Decimal(intermed), roundingMode: RoundingMode.HalfUp)

    if (x_.toFloat64() == 0.0) {
        return 0.0
    }

    return (x_ / pow).toFloat64()
}

public func roundEven(x: Float64, prec: Int64): Float64 {
    if (x == 0.0) {
        return 0.0
    }
    if (prec >= 0 && x == trunc(x)) {
        return x
    }
    let pow = Decimal(10).powWithPrecision(prec, 0)
    let intermed = x * pow.toFloat64()
    if (intermed.isInf()) {
        return x
    }
    if (intermed.isNaN()) {
        return intermed
    }

    let x_ = round(Decimal(intermed), roundingMode: RoundingMode.HalfEven)

    if (x_.toFloat64() == 0.0) {
        return 0.0
    }

    return (x_ / pow).toFloat64()
}

public func same(a: Float64, b: Float64): Bool {
    return a == b || (a.isNaN() && b.isNaN())
}
