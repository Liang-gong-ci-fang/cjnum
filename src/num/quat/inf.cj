package cjnum.num.quat

import std.math.*

// IsInf returns true if any of real(q), imag(q), jmag(q), or kmag(q) is an infinity.
public func nqIsInf(q: Number): Bool {
    return q.real.isInf() || q.imag.isInf() || q.jmag.isInf() || q.kmag.isInf()
}

// Inf returns a quaternion infinity, quaternion(+Inf, +Inf, +Inf, +Inf).
public func nqInf(): Number {
    let inf = Float64.Inf
    return Number(inf, inf, inf, inf)
}
