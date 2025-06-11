package cjnum.num.quat

import cjnum.util.*
import std.math.*

// Sin returns the sine of q.
public func nqSin(q: NqNumber): NqNumber {
    let (w, uv) = nqSplit(q)
    if (uv.nqEqual(zero)) {
        return nqLift(sin(w))
    }
    let v = nqAbs(uv)
    let (s, c) = sincos(w)
    let (sh, ch) = nqSinhcosh(v)
    return nqJoin(s * ch, nqScale(c * sh / v, uv))
}

// Sinh returns the hyperbolic sine of q.
public func nqSinh(q: NqNumber): NqNumber {
    let (w, uv) = nqSplit(q)
    if (uv.nqEqual(zero)) {
        return nqLift(sinh(w))
    }
    let v = nqAbs(uv)
    let (s, c) = sincos(v)
    let (sh, ch) = sinhcosh(w)
    return nqJoin(c * sh, scale(s * ch / v, uv))
}

// Cos returns the cosine of q.
public func nqCos(q: NqNumber): NqNumber {
    let (w, uv) = nqSplit(q)
    if (uv.nqEqual(zero)) {
        return nqLift(cos(w))
    }
    let v = nqAbs(uv)
    let (s, c) = sincos(w)
    let (sh, ch) = sinhcosh(v)
    return nqJoin(c * ch, nqScale(-s * sh / v, uv))
}

// Cosh returns the hyperbolic cosine of q.
public func nqCosh(q: NqNumber): NqNumber {
    let (w, uv) = nqSplit(q)
    if (uv.nqEqual(zero)) {
        return nqLift(cosh(w))
    }
    let v = nqAbs(uv)
    let (s, c) = sincos(v)
    let (sh, ch) = sinhcosh(w)
    return nqJoin(c * ch, scale(s * sh / v, uv))
}

// Tan returns the tangent of q.
public func nqTan(q: NqNumber): NqNumber {
    let d = nqCos(q)
    if (d.nqEqual(zero)) {
        return nqInf()
    }
    return nqMul(nqSin(q), nqInv(d))
}

// Tanh returns the hyperbolic tangent of q.
public func nqTanh(q: NqNumber): NqNumber {
    if (q.real.isInf() && q.real > 0.0) {
        var r = NqNumber(1.0, 0.0, 0.0, 0.0)
        // Change signs dependent on imaginary parts.
        r.imag *= sin(2.0 * q.imag)
        r.jmag *= sin(2.0 * q.jmag)
        r.kmag *= sin(2.0 * q.kmag)
        return r
    }
    let d = nqCosh(q)
    if (d.nqEqual(zero)) {
        return nqInf()
    }
    return nqMul(nqSinh(q), nqInv(d))
}

// Asin returns the inverse sine of q.
public func nqAsin(q: NqNumber): NqNumber {
    let (_, uv) = nqSplit(q)
    if (uv.nqEqual(zero)) {
        return nqLift(asin(q.real))
    }
    let u = nqUnit(uv)
    return nqMul(nqScale(-1.0, u), nqLog(nqAdd(nqMul(u, q), nqSqrt(nqSub(NqNumber(1.0, 0.0, 0.0, 0.0), nqMul(q, q))))))
}

// Asinh returns the inverse hyperbolic sine of q.
public func nqAsinh(q: NqNumber): NqNumber {
    return nqLog(nqAdd(q, nqSqrt(nqAdd(NqNumber(1.0, 0.0, 0.0, 0.0), nqMul(q, q)))))
}

// Acos returns the inverse cosine of q.
public func nqAcos(q: NqNumber): NqNumber {
    let (w, uv) = nqSplit(nqAsin(q))
    return nqJoin(Float64.PI / 2.0 - w, nqScale(-1.0, uv))
}

// Acosh returns the inverse hyperbolic cosine of q.
public func nqAcosh(q: NqNumber): NqNumber {
    var w = nqAcos(q)
    let (_, uv) = nqSplit(w)
    if (uv.nqEqual(zero)) {
        return w
    }
    w = nqMul(w, nqUnit(uv))
    if (w.real < 0.0) {
        w = nqScale(-1.0, w)
    }
    return w
}

// Atan returns the inverse tangent of q.
public func nqAtan(q: NqNumber): NqNumber {
    let (w, uv) = nqSplit(q)
    if (uv.nqEqual(zero)) {
        return nqLift(atan(w))
    }
    let u = nqUnit(uv)
    return nqMul(nqMul(nqLift(0.5), u), nqLog(nqMul(nqAdd(u, q), nqInv(nqSub(u, q)))))
}

// Atanh returns the inverse hyperbolic tangent of q.
public func nqAtanh(q: NqNumber): NqNumber {
    let (w, uv) = nqSplit(q)
    if (uv.nqEqual(zero)) {
        return nqLift(nanAtanh(w))
    }
    let u = nqUnit(uv)
    return nqMul(nqScale(-1.0, u), nqAtan(nqMul(u, q)))
}

// calculate sinh and cosh
public func nqSinhcosh(x: Float64): (Float64, Float64) {
    if (abs(x) <= 0.5) {
        return (sinh(x), cosh(x))
    }
    var e = exp(x)
    let ei = 0.5 / e
    e *= 0.5
    return (e - ei, e + ei)
}

// scale returns q scaled by f, except that inf×0 is 0.
private func scale(f: Float64, q: NqNumber): NqNumber {
    if (f == 0.0) {
        return NqNumber(0.0, 0.0, 0.0, 0.0)
    }
    if (q.real != 0.0) {
        q.real *= f
    }
    if (q.imag != 0.0) {
        q.imag *= f
    }
    if (q.jmag != 0.0) {
        q.jmag *= f
    }
    if (q.kmag != 0.0) {
        q.kmag *= f
    }
    return q
}
