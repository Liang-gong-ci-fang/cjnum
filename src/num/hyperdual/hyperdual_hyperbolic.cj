package cjnum.num.hyperdual

import cjnum.util.*
import std.math.*

// Sinh returns the hyperbolic sine of d.
//
// Special cases are:
//
//	Sinh(±0) = (±0+Nϵ₁+Nϵ₂±0ϵ₁ϵ₂)
//	Sinh(±Inf) = ±Inf
//	Sinh(NaN) = NaN
public func nhSinh(d: NhNumber): NhNumber {
    if (d.real == 0.0) {
        return NhNumber(
            d.real,
            d.e1mag,
            d.e1mag,
            d.real,
        )
    }
    if (d.real.isInf()) {
        return NhNumber(
            d.real,
            Float64.Inf,
            Float64.Inf,
            d.real,
        )
    }
    let fn = nanSinh(d.real)
    let deriv = nanCosh(d.real)
    return NhNumber(
        fn,
        deriv * d.e1mag,
        deriv * d.e2mag,
        deriv * d.e1E2mag + fn * d.e1mag * d.e2mag,
    )
}

// Cosh returns the hyperbolic cosine of d.
//
// Special cases are:
//
//	Cosh(±0) = 1
//	Cosh(±Inf) = +Inf
//	Cosh(NaN) = NaN
public func nhCosh(d: NhNumber): NhNumber {
    if (d.real.isInf()) {
        return NhNumber(
            Float64.Inf,
            d.real,
            d.real,
            Float64.Inf,
        )
    }
    let fn = nanCosh(d.real)
    let deriv = nanSinh(d.real)
    return NhNumber(
        fn,
        deriv * d.e1mag,
        deriv * d.e2mag,
        deriv * d.e1E2mag + fn * d.e1mag * d.e2mag,
    )
}

// Tanh returns the hyperbolic tangent of d.
//
// Special cases are:
//
//	Tanh(±0) = (±0+Nϵ₁+Nϵ₂∓0ϵ₁ϵ₂)
//	Tanh(±Inf) = (±1+0ϵ₁+0ϵ₂∓0ϵ₁ϵ₂)
//	Tanh(NaN) = NaN
public func nhTanh(d: NhNumber): NhNumber {
    if (d.real == 0.0) {
        return NhNumber(
            d.real,
            d.e1mag,
            d.e2mag,
            -d.real,
        )
    } else if (d.real == Float64.Inf) {
        return NhNumber(
            1.0,
            0.0,
            0.0,
            negZero,
        )
    } else if (d.real == -Float64.Inf) {
        return NhNumber(
            -1.0,
            0.0,
            0.0,
            0.0,
        )
    }
    let fn = nanTanh(d.real)
    let deriv = 1.0 - fn * fn
    return NhNumber(
        fn,
        deriv * d.e1mag,
        deriv * d.e2mag,
        deriv * d.e1E2mag - d.e1mag * d.e2mag * (2.0 * fn * deriv),
    )
}

// Asinh returns the inverse hyperbolic sine of d.
//
// Special cases are:
//
//	Asinh(±0) = (±0+Nϵ₁+Nϵ₂∓0ϵ₁ϵ₂)
//	Asinh(±Inf) = ±Inf
//	Asinh(NaN) = NaN
public func nhAsinh(d: NhNumber): NhNumber {
    if (d.real == 0.0) {
        return NhNumber(
            d.real,
            d.e1mag,
            d.e2mag,
            -d.real,
        )
    }
    let fn = asinh(d.real)
    let deriv1 = d.real * d.real + 1.0
    let deriv = 1.0 / sqrt(deriv1)
    return NhNumber(
        fn,
        deriv * d.e1mag,
        deriv * d.e2mag,
        deriv * d.e1E2mag + d.e1mag * d.e2mag * (-d.real * (deriv / deriv1)),
    )
}

// Acosh returns the inverse hyperbolic cosine of d.
//
// Special cases are:
//
//	Acosh(+Inf) = +Inf
//	Acosh(1) = (0+Infϵ₁+Infϵ₂-Infϵ₁ϵ₂)
//	Acosh(x) = NaN if x < 1
//	Acosh(NaN) = NaN
public func nhAcosh(d: NhNumber): NhNumber {
    if (d.real <= 1.0) {
        if (d.real == 1.0) {
            return NhNumber(
                0.0,
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
    let fn = nanAcosh(d.real)
    let deriv1 = d.real * d.real - 1.0
    let deriv = 1.0 / nanSqrt(deriv1)
    return NhNumber(
        fn,
        deriv * d.e1mag,
        deriv * d.e2mag,
        deriv * d.e1E2mag + d.e1mag * d.e2mag * (-d.real * (deriv / deriv1)),
    )
}

// Atanh returns the inverse hyperbolic tangent of d.
//
// Special cases are:
//
//	Atanh(1) = +Inf
//	Atanh(±0) = (±0+Nϵ₁+Nϵ₂±0ϵ₁ϵ₂)
//	Atanh(-1) = -Inf
//	Atanh(x) = NaN if x < -1 or x > 1
//	Atanh(NaN) = NaN
public func nhAtanh(d: NhNumber): NhNumber {
    if (d.real == 0.0) {
        return NhNumber(
            d.real,
            d.e1mag,
            d.e2mag,
            d.real,
        )
    }
    if (abs(d.real) == 1.0) {
        if (d.real > 0.0) {
            return NhNumber(
                Float64.Inf,
                Float64.NaN,
                Float64.NaN,
                Float64.Inf,
            )
        } else {
            return NhNumber(
                -Float64.Inf,
                Float64.NaN,
                Float64.NaN,
                -Float64.Inf,
            )
        }
    }
    let fn = nanAtanh(d.real)
    let deriv1 = 1.0 - d.real * d.real
    let deriv = 1.0 / deriv1
    return NhNumber(
        fn,
        deriv * d.e1mag,
        deriv * d.e2mag,
        deriv * d.e1E2mag + d.e1mag * d.e2mag * (2.0 * d.real / (deriv1 * deriv1)),
    )
}
