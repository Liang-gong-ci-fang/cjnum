package cjnum.num.dual

import cjnum.util.*
import std.math.*

// PowReal returns x**p, the base-x exponential of p.
//
// Special cases are (in order):
//
//	PowReal(NaN+xϵ, ±0) = 1+NaNϵ for any x
//	PowReal(x, ±0) = 1 for any x
//	PowReal(1+xϵ, y) = 1+xyϵ for any y
//	PowReal(x, 1) = x for any x
//	PowReal(NaN+xϵ, y) = NaN+NaNϵ
//	PowReal(x, NaN) = NaN+NaNϵ
//	PowReal(±0, y) = ±Inf for y an odd integer < 0
//	PowReal(±0, -Inf) = +Inf
//	PowReal(±0, +Inf) = +0
//	PowReal(±0, y) = +Inf for finite y < 0 and not an odd integer
//	PowReal(±0, y) = ±0 for y an odd integer > 0
//	PowReal(±0, y) = +0 for finite y > 0 and not an odd integer
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
//	PowReal(x, y) = NaN+NaNϵ for finite x < 0 and finite non-integer y
public func ndPowReal(d: Number, p: Float64): Number {
    let tol = 1e-15

    var r = d.real
    if (abs(r) < tol) {
        if (r >= 0.0) {
            r = tol
        }
        if (r < 0.0) {
            r = -tol
        }
    }
    let deriv = p * pow(r, p - 1.0)
    return Number(
        pow(d.real, p),
        d.emag * deriv,
    )
}

// Sqrt returns the square root of d.
//
// Special cases are:
//
//	Sqrt(+Inf) = +Inf
//	Sqrt(±0) = (±0+Infϵ)
//	Sqrt(x < 0) = NaN
//	Sqrt(NaN) = NaN
public func ndSqrt(d: Number): Number {
    if (d.real <= 0.0) {
        if (d.real == 0.0) {
            return Number(
                d.real,
                Float64.Inf,
            )
        }
        return Number(
            Float64.NaN,
            Float64.NaN,
        )
    }
    return ndPowReal(d, 0.5)
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
public func ndExp(d: Number): Number {
    let fnDeriv = exp(d.real)
    return Number(
        fnDeriv,
        fnDeriv * d.emag,
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
public func ndLog(d: Number): Number {
    if (d.real == 0.0) {
        return Number(
            log(d.real),
            copysign(Float64.Inf, d.real)
        )
    } else if (d.real == Float64.Inf) {
        return Number(
            log(d.real),
            0.0,
        )
    }
    if (d.real < 0.0) {
        return Number(
            Float64.NaN,
            Float64.NaN,
        )
    }
    return Number(
        log(d.real),
        d.emag / d.real,
    )
}

// Sin returns the sine of d.
//
// Special cases are:
//
//	Sin(±0) = (±0+Nϵ)
//	Sin(±Inf) = NaN
//	Sin(NaN) = NaN
public func ndSin(d: Number): Number {
    if (d.real == 0.0) {
        return Number(
            d.real,
            d.emag,
        )
    }
    let fn = sin(d.real)
    let deriv = cos(d.real)
    return Number(
        fn,
        deriv * d.emag,
    )
}

// Cos returns the cosine of d.
//
// Special cases are:
//
//	Cos(±Inf) = NaN
//	Cos(NaN) = NaN
public func ndCos(d: Number): Number {
    let fn = cos(d.real)
    let deriv = -sin(d.real)
    return Number(
        fn,
        deriv * d.emag,
    )
}

// Tan returns the tangent of d.
//
// Special cases are:
//
//	Tan(±0) = (±0+Nϵ)
//	Tan(±Inf) = NaN
//	Tan(NaN) = NaN
public func ndTan(d: Number): Number {
    if (d.real == 0.0) {
        return Number(
            d.real,
            d.emag,
        )
    }
    let fn = tan(d.real)
    let deriv = 1.0 + fn * fn
    return Number(
        fn,
        deriv * d.emag,
    )
}

// Asin returns the inverse sine of d.
//
// Special cases are:
//
//	Asin(±0) = (±0+Nϵ)
//	Asin(±1) = (±Inf+Infϵ)
//	Asin(x) = NaN if x < -1 or x > 1
public func ndAsin(d: Number): Number {
    let m = abs(d.real);
    if (d.real == 0.0) {
        return Number(
            d.real,
            d.emag,
        )
    } else if (m >= 1.0) {
        if (m == 1.0) {
            return Number(
                asin(d.real),
                Float64.Inf,
            )
        }
        return Number(
            Float64.NaN,
            Float64.NaN,
        )
    }
    let fn = asin(d.real)
    let deriv = 1.0 / sqrt(1.0 - d.real * d.real)
    return Number(
        fn,
        deriv * d.emag,
    )
}

// Acos returns the inverse cosine of d.
//
// Special cases are:
//
//	Acos(-1) = (Pi-Infϵ)
//	Acos(1) = (0-Infϵ)
//	Acos(x) = NaN if x < -1 or x > 1
public func ndAcos(d: Number): Number {
    let m = abs(d.real)
    if (m >= 1.0) {
        if (m == 1.0) {
            return Number(
                acos(d.real),
                -Float64.Inf,
            )
        }
        return Number(
            Float64.NaN,
            Float64.NaN,
        )
    }
    let fn = acos(d.real)
    let deriv = -1.0 / sqrt(1.0 - d.real * d.real)
    return Number(
        fn,
        deriv * d.emag,
    )
}

// Atan returns the inverse tangent of d.
//
// Special cases are:
//
//	Atan(±0) = (±0+Nϵ)
//	Atan(±Inf) = (±Pi/2+0ϵ)
public func ndAtan(d: Number): Number {
    if (d.real == 0.0) {
        return Number(
            d.real,
            d.emag,
        )
    }
    let fn = atan(d.real)
    let deriv = 1.0 / (1.0 + d.real * d.real)
    return Number(
        fn,
        deriv * d.emag,
    )
}
