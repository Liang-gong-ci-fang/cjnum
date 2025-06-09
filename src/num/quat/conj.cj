package cjnum.num.quat

import std.math.*

// Conj returns the quaternion conjugate of q.
public func nqConj(q: Number): Number {
    return Number(q.real, -q.imag, -q.jmag, -q.kmag)
}

// Inv returns the quaternion inverse of q.
public func nqInv(q: Number): Number {
    if (nqIsInf(q)) {
        return zero
    }
    let a = nqAbs(q)
    return nqScale(1.0 / (a * a), nqConj(q))
}
