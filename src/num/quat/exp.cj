package cjnum.num.quat

import cjnum.util.*
import std.math.*

// Exp returns e**q, the base-e exponential of q.
public func nqExp(q: Number): Number {
    var (w, uv) = nqSplit(q)
    if (uv.nqEqual(zero)) {
        return nqLift(exp(w))
    }
    let v = nqAbs(uv)
    let e = exp(w)
    let (s, c) = sincos(v)
    return nqJoin(e * c, nqScale(e * s / v, uv))
}

// Log returns the natural logarithm of q.
public func nqLog(q: Number): Number {
    var (w, uv) = nqSplit(q)
    if (uv.nqEqual(zero)) {
        return nqLift(log(w))
    }
    let v = nqAbs(uv)
    return nqJoin(log(nqAbs(q)), nqScale(atan2(v, w) / v, uv))
}

// Pow return q**r, the base-q exponential of r.
// For generalized compatibility with math.Pow:
//
//	Pow(0, ±0) returns 1+0i+0j+0k
//	Pow(0, c) for real(c)<0 returns Inf+0i+0j+0k if imag(c), jmag(c), kmag(c) are zero,
//	    otherwise Inf+Inf i+Inf j+Inf k.
public func nqPow(q: Number, r: Number): Number {
    if (q.nqEqual(zero)) {
        var (w, uv) = nqSplit(r)
        if (w == 0.0) {
            return Number(1.0, 0.0, 0.0, 0.0)
        } else if (w < 0.0) {
            if (uv.nqEqual(zero)) {
                return Number(Float64.Inf, 0.0, 0.0, 0.0)
            }
            return nqInf()
        } else if (w > 0.0) {
            return zero
        }
    }
    return nqExp(nqMul(nqLog(q), r))
}

// PowReal return q**r, the base-q exponential of r.
// For generalized compatibility with math.Pow:
//
//	PowReal(0, ±0) returns 1+0i+0j+0k
//	PowReal(0, c) for c<0 returns Inf+0i+0j+0k.
public func nqPowReal(q: Number, r: Float64): Number {
    if (q.nqEqual(zero)) {
        if (r == 0.0) {
            return Number(1.0, 0.0, 0.0, 0.0)
        } else if (r < 0.0) {
            return nqInf()
        } else if (r > 0.0) {
            return zero
        }
    }
    return nqExp(nqScale(r, nqLog(q)))
}

// Sqrt returns the square root of q.
public func nqSqrt(q: Number): Number {
    if (q.nqEqual(zero)) {
        return zero
    }
    return nqPowReal(q, 0.5)
}
