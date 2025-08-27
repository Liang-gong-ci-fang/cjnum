package cjnum.num.quat

// IsInf returns true if any of real(q), imag(q), jmag(q), or kmag(q) is an infinity.
public func nqIsInf(q: NqNumber): Bool {
    return q.real.isInf() || q.imag.isInf() || q.jmag.isInf() || q.kmag.isInf()
}

// Inf returns a quaternion infinity, quaternion(+Inf, +Inf, +Inf, +Inf).
public func nqInf(): NqNumber {
    let inf = Float64.Inf
    return NqNumber(inf, inf, inf, inf)
}
