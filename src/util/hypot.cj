package cjnum.util

import std.math.*

/*
Hypot -- sqrt(p*p + q*q), but overflows only if the result does.
 */

// Hypot returns [Sqrt](p*p + q*q), taking care to avoid
// unnecessary overflow and underflow.
//
// Special cases are:
//
//	Hypot(±Inf, q) = +Inf
//	Hypot(p, ±Inf) = +Inf
//	Hypot(NaN, q) = NaN
//	Hypot(p, NaN) = NaN
public func hypot(p: Float64, q: Float64): Float64 {
    // TODO if haveArchHypot { 
    return _hypot(p, q)
}

private func _hypot(p_: Float64, q_: Float64): Float64 {
    var p = p_
    var q = q_
    p = abs(p)
    q = abs(q)
    // special cases
    if ((p.isInf()) || (q.isInf())) {
        return Float64.Inf
    } else if (p.isNaN() || q.isNaN()) {
        return Float64.NaN
    }
    if (p < q) {
        var temp = p
        p = q
        q = temp
    }
    if (p == 0.0) {
        return 0.0
    }
    q = q / p
    return p * sqrt(1.0 + q * q)
}
