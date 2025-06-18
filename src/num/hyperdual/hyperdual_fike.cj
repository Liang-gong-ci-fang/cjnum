package cjnum.num.hyperdual

import cjnum.util.*
import std.math.*

// PowReal returns x**p, the base-x exponential of p.
//
// Special cases are (in order):
//
//	PowReal(NaN+xϵ₁+yϵ₂, ±0) = 1+NaNϵ₁+NaNϵ₂+NaNϵ₁ϵ₂ for any x and y
//	PowReal(x, ±0) = 1 for any x
//	PowReal(1+xϵ₁+yϵ₂, z) = 1+xzϵ₁+yzϵ₂+2xyzϵ₁ϵ₂ for any z
//	PowReal(NaN+xϵ₁+yϵ₂, 1) = NaN+xϵ₁+yϵ₂+NaNϵ₁ϵ₂ for any x
//	PowReal(x, 1) = x for any x
//	PowReal(NaN+xϵ₁+xϵ₂, y) = NaN+NaNϵ₁+NaNϵ₂+NaNϵ₁ϵ₂
//	PowReal(x, NaN) = NaN+NaNϵ₁+NaNϵ₂+NaNϵ₁ϵ₂
//	PowReal(±0, y) = ±Inf for y an odd integer < 0
//	PowReal(±0, -Inf) = +Inf
//	PowReal(±0, +Inf) = +0
//	PowReal(±0, y) = +Inf for finite y < 0 and not an odd integer
//	PowReal(±0, y) = ±0 for y an odd integer > 0
//	PowReal(±0, y) = +0 for finite y > 0 and not an odd integer
//	PowReal(-1, ±Inf) = 1
//	PowReal(x+0ϵ₁+0ϵ₂, +Inf) = +Inf+NaNϵ₁+NaNϵ₂+NaNϵ₁ϵ₂ for |x| > 1
//	PowReal(x+xϵ₁+yϵ₂, +Inf) = +Inf+Infϵ₁+Infϵ₂+NaNϵ₁ϵ₂ for |x| > 1
//	PowReal(x, -Inf) = +0+NaNϵ₁+NaNϵ₂+NaNϵ₁ϵ₂ for |x| > 1
//	PowReal(x+yϵ₁+zϵ₂, +Inf) = +0+NaNϵ₁+NaNϵ₂+NaNϵ₁ϵ₂ for |x| < 1
//	PowReal(x+0ϵ₁+0ϵ₂, -Inf) = +Inf+NaNϵ₁+NaNϵ₂+NaNϵ₁ϵ₂ for |x| < 1
//	PowReal(x, -Inf) = +Inf-Infϵ₁-Infϵ₂+NaNϵ₁ϵ₂ for |x| < 1
//	PowReal(+Inf, y) = +Inf for y > 0
//	PowReal(+Inf, y) = +0 for y < 0
//	PowReal(-Inf, y) = Pow(-0, -y)
//	PowReal(x, y) = NaN+NaNϵ₁+NaNϵ₂+NaNϵ₁ϵ₂ for finite x < 0 and finite non-integer y
public func nhPowReal(d: NhNumber, p: Float64): NhNumber {
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
    return NhNumber(
        pow(d.real, p),
        d.e1mag * deriv,
        d.e2mag * deriv,
        d.e1E2mag * deriv + p * (p - 1.0) * d.e1mag * d.e2mag * pow(r, (p - 2.0)),
    )
}

// Pow returns x**p, the base-x exponential of p.
public func nhPow(d: NhNumber, p: NhNumber): NhNumber {
    return nhExp(nhMul(p, nhLog(d)))
}

// Sqrt returns the square root of d.
//
// Special cases are:
//
//	Sqrt(+Inf) = +Inf
//	Sqrt(±0) = (±0+Infϵ₁+Infϵ₂-Infϵ₁ϵ₂)
//	Sqrt(x < 0) = NaN
//	Sqrt(NaN) = NaN
public func nhSqrt(d: NhNumber): NhNumber {
    if (d.real <= 0.0) {
        if (d.real == 0.0) {
            return NhNumber(
                d.real,
                Float64.Inf,
                Float64.Inf,
                -Float64.Inf,
            )
        }
        return NhNumber(
            Float64.NaN,
            Float64.NaN,
            Float64.NaN,
            Float64.NaN,
        )
    }
    return nhPowReal(d, 0.5)
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
public func nhExp(d: NhNumber): NhNumber {
    let exp = exp(d.real) // exp is also the derivative.
    return NhNumber(
        exp,
        exp * d.e1mag,
        exp * d.e2mag,
        exp * (d.e1E2mag + d.e1mag * d.e2mag),
    )
}

// Log returns the natural logarithm of d.
//
// Special cases are:
//
//	Log(+Inf) = (+Inf+0ϵ₁+0ϵ₂-0ϵ₁ϵ₂)
//	Log(0) = (-Inf±Infϵ₁±Infϵ₂-Infϵ₁ϵ₂)
//	Log(x < 0) = NaN
//	Log(NaN) = NaN
public func nhLog(d: NhNumber): NhNumber {
    if (d.real == 0.0) {
        return NhNumber(
            log(d.real),
            copysign(Float64.Inf, d.real),
            copysign(Float64.Inf, d.real),
            -Float64.Inf,
        )
    } else if (d.real == Float64.Inf) {
        return NhNumber(
            log(d.real),
            0.0,
            0.0,
            negZero,
        )
    } else if (d.real < 0.0) {
        return NhNumber(
            Float64.NaN,
            Float64.NaN,
            Float64.NaN,
            Float64.NaN,
        )
    }
    let deriv1 = d.e1mag / d.real
    let deriv2 = d.e2mag / d.real
    return NhNumber(
        log(d.real),
        deriv1,
        deriv2,
        d.e1E2mag / d.real - (deriv1 * deriv2),
    )
}

// Sin returns the sine of d.
//
// Special cases are:
//
//	Sin(±0) = (±0+Nϵ₁+Nϵ₂∓0ϵ₁ϵ₂)
//	Sin(±Inf) = NaN
//	Sin(NaN) = NaN
public func nhSin(d: NhNumber): NhNumber {
    if (d.real == 0.0) {
        return NhNumber(
            d.real,
            d.e1mag,
            d.e2mag,
            -d.real,
        )
    }
    let fn = sin(d.real)
    let deriv = cos(d.real)
    return NhNumber(
        fn,
        deriv * d.e1mag,
        deriv * d.e2mag,
        deriv * d.e1E2mag - fn * d.e1mag * d.e2mag,
    )
}

// Cos returns the cosine of d.
//
// Special cases are:
//
//	Cos(±Inf) = NaN
//	Cos(NaN) = NaN
public func nhCos(d: NhNumber): NhNumber {
    let fn = cos(d.real)
    let deriv = -sin(d.real)
    return NhNumber(
        fn,
        deriv * d.e1mag,
        deriv * d.e2mag,
        deriv * d.e1E2mag - fn * d.e1mag * d.e2mag,
    )
}

// Tan returns the tangent of d.
//
// Special cases are:
//
//	Tan(±0) = (±0+Nϵ₁+Nϵ₂±0ϵ₁ϵ₂)
//	Tan(±Inf) = NaN
//	Tan(NaN) = NaN
public func nhTan(d: NhNumber): NhNumber {
    if (d.real == 0.0) {
        return NhNumber(
            d.real,
            d.e1mag,
            d.e2mag,
            d.real,
        )
    }
    let fn = tan(d.real)
    let deriv = 1.0 + fn * fn
    return NhNumber(
        fn,
        deriv * d.e1mag,
        deriv * d.e2mag,
        deriv * d.e1E2mag + d.e1mag * d.e2mag * (2.0 * fn * deriv),
    )
}

// Asin returns the inverse sine of d.
//
// Special cases are:
//
//	Asin(±0) = (±0+Nϵ₁+Nϵ₂±0ϵ₁ϵ₂)
//	Asin(±1) = (±Inf+Infϵ₁+Infϵ₂±Infϵ₁ϵ₂)
//	Asin(x) = NaN if x < -1 or x > 1
public func nhAsin(d: NhNumber): NhNumber {
    let m = abs(d.real);
    if (d.real == 0.0) {
        return NhNumber(
            d.real,
            d.e1mag,
            d.e2mag,
            d.real,
        )
    } else if (m >= 1.0) {
        if (m == 1.0) {
            return NhNumber(
                asin(d.real),
                Float64.Inf,
                Float64.Inf,
                copysign(Float64.Inf, d.real),
            )
        }
        return NhNumber(
            Float64.NaN,
            Float64.NaN,
            Float64.NaN,
            Float64.NaN,
        )
    }
    let fn = nanAsin(d.real)
    let deriv1 = 1.0 - d.real * d.real
    let deriv = 1.0 / sqrt(deriv1)
    return NhNumber(
        fn,
        deriv * d.e1mag,
        deriv * d.e2mag,
        deriv * d.e1E2mag + d.e1mag * d.e2mag * (d.real * pow(deriv1, -1.5)),
    )
}

// Acos returns the inverse cosine of d.
//
// Special cases are:
//
//	Acos(-1) = (Pi-Infϵ₁-Infϵ₂+Infϵ₁ϵ₂)
//	Acos(1) = (0-Infϵ₁-Infϵ₂-Infϵ₁ϵ₂)
//	Acos(x) = NaN if x < -1 or x > 1
public func nhAcos(d: NhNumber): NhNumber {
    let m = abs(d.real);
    if (m >= 1.0) {
        if (m == 1.0) {
            return NhNumber(
                acos(d.real),
                -Float64.Inf,
                -Float64.Inf,
                copysign(Float64.Inf, -d.real),
            )
        }
        return NhNumber(
            Float64.NaN,
            Float64.NaN,
            Float64.NaN,
            Float64.NaN,
        )
    }
    let fn = nanAcos(d.real)
    let deriv1 = 1.0 - d.real * d.real
    let deriv = -1.0 / sqrt(deriv1)
    return NhNumber(
        fn,
        deriv * d.e1mag,
        deriv * d.e2mag,
        deriv * d.e1E2mag + d.e1mag * d.e2mag * (-d.real * pow(deriv1, -1.5)),
    )
}

// Atan returns the inverse tangent of d.
//
// Special cases are:
//
//	Atan(±0) = (±0+Nϵ₁+Nϵ₂∓0ϵ₁ϵ₂)
//	Atan(±Inf) = (±Pi/2+0ϵ₁+0ϵ₂∓0ϵ₁ϵ₂)
public func nhAtan(d: NhNumber): NhNumber {
    if (d.real == 0.0) {
        return NhNumber(
            d.real,
            d.e1mag,
            d.e2mag,
            -d.real,
        )
    }
    let fn = atan(d.real)
    let deriv1 = 1.0 + d.real * d.real
    let deriv = 1.0 / deriv1
    return NhNumber(
        fn,
        deriv * d.e1mag,
        deriv * d.e2mag,
        deriv * d.e1E2mag + d.e1mag * d.e2mag * (-2.0 * d.real / (deriv1 * deriv1)),
    )
}
