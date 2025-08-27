package cjnum.num.quat


// IsNaN returns true if any of real(q), imag(q), jmag(q), or kmag(q) is NaN
// and none are an infinity.
public func nqIsNaN(q: NqNumber): Bool {
    if (q.real.isInf() || q.imag.isInf() || q.jmag.isInf() || q.kmag.isInf()) {
        return false
    }
    return q.real.isNaN() || q.imag.isNaN() || q.jmag.isNaN() || q.kmag.isNaN()
}

// NaN returns a quaternion “not-a-number” value.
public func nqNaN(): NqNumber {
    let nan = Float64.NaN
    return NqNumber(nan, nan, nan, nan)
}
