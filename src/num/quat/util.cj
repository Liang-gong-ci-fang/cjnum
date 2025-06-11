package cjnum.num.quat

import std.math.*

public func nqLift(v: Float64): NqNumber {
    return NqNumber(v, 0.0, 0.0, 0.0)
}

public func nqSplit(q: NqNumber): (Float64, NqNumber) {
    return (q.real, NqNumber(0.0, q.imag, q.jmag, q.kmag))
}

public func nqJoin(w: Float64, uv: NqNumber): NqNumber {
    return NqNumber(w, uv.imag, uv.jmag, uv.kmag)
}

public func nqUnit(q: NqNumber): NqNumber {
    return nqScale(1.0 / nqAbs(q), q)
}
