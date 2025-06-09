package cjnum.num.quat

import std.math.*

public func nqLift(v: Float64): Number {
    return Number(v, 0.0, 0.0, 0.0)
}

public func nqSplit(q: Number): (Float64, Number) {
    return (q.real, Number(0.0, q.imag, q.jmag, q.kmag))
}

public func nqJoin(w: Float64, uv: Number): Number {
    return Number(w, uv.imag, uv.jmag, uv.kmag)
}

public func nqUnit(q: Number): Number {
    return nqScale(1.0 / nqAbs(q), q)
}
