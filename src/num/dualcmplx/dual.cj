package cjnum.num.dualcmplx

import cjnum.complex.*
import cjnum.util.*
import std.math.*

// NdcNumber is a float64 precision anti-commutative dual complex number.
public class NdcNumber {
    let real: Complex128
    let dual: Complex128
    public init(real: Complex128, dual: Complex128) {
        this.real = real
        this.dual = dual
    }
}

// Abs returns the absolute value of d.
public func ndcAbs(d: NdcNumber): Float64 {
    return c128Abs(d.real)
}

// PowReal returns d**p, the base-d exponential of p.
//
// Special cases are (in order):
//
//	PowReal(NaN+xϵ, ±0) = 1+NaNϵ for any x
//	Pow(0+xϵ, y) = 0+Infϵ for all y < 1.
//	Pow(0+xϵ, y) = 0 for all y > 1.
//	PowReal(x, ±0) = 1 for any x
//	PowReal(1+xϵ, y) = 1+xyϵ for any y
//	Pow(Inf, y) = +Inf+NaNϵ for y > 0
//	Pow(Inf, y) = +0+NaNϵ for y < 0
//	PowReal(x, 1) = x for any x
//	PowReal(NaN+xϵ, y) = NaN+NaNϵ
//	PowReal(x, NaN) = NaN+NaNϵ
//	PowReal(-1, ±Inf) = 1
//	PowReal(x+0ϵ, +Inf) = +Inf+NaNϵ for |x| > 1
//	PowReal(x+yϵ, +Inf) = +Inf for |x| > 1
//	PowReal(x, -Inf) = +0+NaNϵ for |x| > 1
//	PowReal(x, +Inf) = +0+NaNϵ for |x| < 1
//	PowReal(x+0ϵ, -Inf) = +Inf+NaNϵ for |x| < 1
//	PowReal(x, -Inf) = +Inf-Infϵ for |x| < 1
//	PowReal(+Inf, y) = +Inf for y > 0
//	PowReal(+Inf, y) = +0 for y < 0
//	PowReal(-Inf, y) = Pow(-0, -y)
public func ndcPowReal(d_: NdcNumber, p: Float64): NdcNumber {
    var d = d_
    if (p == 0.0) {
        if (d.real.isNaN()) {
            return NdcNumber(Complex128(1.0, 0.0), Complex128.NaN)
        } else if ((d.real.real == 0.0 && d.real.imag == 0.0) || d.real.isInf()) {
            return NdcNumber(Complex128(1.0, 0.0), Complex128(0.0, 0.0))
        }
    } else if (p == 1.0) {
        if (d.real.isInf()) {
            return NdcNumber(d.real, Complex128.NaN)
        }
        return d
    } else if (p.isInf() && p > 0.0) {
        if (d.real.real == -1.0 && d.real.imag == 0.0) {
            return NdcNumber(Complex128(1.0, 0.0), Complex128.NaN)
        }
        if (ndcAbs(d) > 1.0) {
            if (d.dual.real == 0.0 && d.dual.imag == 0.0) {
                return NdcNumber(Complex128.Inf, Complex128.NaN)
            }
            return NdcNumber(Complex128.Inf, Complex128.Inf)
        }
        return NdcNumber(Complex128(0.0, 0.0), Complex128.NaN)
    } else if (p.isInf() && p < 0.0) {
        if (d.real.real == -1.0 && d.real.imag == 0.0) {
            return NdcNumber(Complex128(1.0, 0.0), Complex128.NaN)
        }
        if (ndcAbs(d) > 1.0) {
            return NdcNumber(Complex128(0.0, 0.0), Complex128.NaN)
        }
        if (d.dual == Complex128(0.0, 0.0)) {
            return NdcNumber(Complex128.Inf, Complex128.NaN)
        }
        return NdcNumber(Complex128.Inf, Complex128.Inf)
    } else if (p.isNaN()) {
        return NdcNumber(Complex128.NaN, Complex128.NaN)
    } else if (d.real == Complex128(0.0, 0.0)) {
        if (p < 1.0) {
            return NdcNumber(d.real, Complex128.Inf)
        }
        return NdcNumber(d.real, Complex128(0.0, 0.0))
    } else if (d.real.isInf()) {
        if (p < 0.0) {
            return NdcNumber(Complex128(0.0, 0.0), Complex128.NaN)
        }
        return NdcNumber(Complex128.Inf, Complex128.NaN)
    }
    return ndcPow(d, NdcNumber(Complex128(p, 0.0), Complex128(0.0, 0.0)))
}

// Pow returns d**p, the base-d exponential of p.
public func ndcPow(d: NdcNumber, p: NdcNumber): NdcNumber {
    return ndcExp(ndcMul(p, ndcLog(d)))
}

// TODO func Sqrt(d NdcNumber) NdcNumber

// Mul returns the dual product of x and y, x×y.
public func ndcMul(x: NdcNumber, y: NdcNumber): NdcNumber {
    return NdcNumber(
        x.real * y.real,
        x.real * y.dual + x.dual * y.real.conj(),
    )
}

// Inv returns the dual inverse of d.
public func ndcInv(d: NdcNumber): NdcNumber {
    return NdcNumber(
        Complex128(1.0, 0.0) / d.real,
        Complex128(-1.0, 0.0) * d.dual / (d.real * d.real.conj()),
    )
}

// Exp returns e**q, the base-e exponential of d.
//
// Special cases are:
//
//	Exp(+Inf) = +Inf
//	Exp(NaN) = NaN
//
// Very large values overflow to 0 or +Inf.
// Very small values underflow to 1.
public func ndcExp(d: NdcNumber): NdcNumber {
    let fn = c128Exp(d.real)
    if (d.real.imag == 0.0) {
        return NdcNumber(fn, fn * d.dual)
    }
    let conj = d.real.conj()
    return NdcNumber(
        fn,
        ((fn - c128Exp(conj)) / (d.real - conj)) * d.dual,
    )
}

// Log returns the natural logarithm of d.
//
// Special cases are:
//
//	Log(+Inf) = (+Inf+0ϵ)
//	Log(0) = (-Inf±Infϵ)
//	Log(x < 0) = NaN
//	Log(NaN) = NaN
public func ndcLog(d: NdcNumber): NdcNumber {
    let fn = c128Log(d.real) //这里把real部分由-0.0 -0.0 转为0.0 0.0
    if (d.real.real == 0.0 && d.real.imag == 0.0) {
        return NdcNumber(
            fn,
            Complex128(copysign(Float64.Inf, d.real.real), Float64.NaN)
        )
    } else if (d.real.imag == 0.0) {
        return NdcNumber(
            fn,
            d.dual / d.real
        )
    } else if (d.real.isInf()) {
        return NdcNumber(
            fn,
            Complex128(0.0, 0.0)
        )
    }
    let conj = d.real.conj()
    return NdcNumber(
        fn,
        ((fn - c128Log(conj)) / (d.real - conj)) * d.dual
    )
}
