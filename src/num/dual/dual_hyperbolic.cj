package cjnum.num.dual

import cjnum.util.*
import std.math.*

// Sinh returns the hyperbolic sine of d.
//
// Special cases are:
//
//	Sinh(±0) = (±0+Nϵ)
//	Sinh(±Inf) = ±Inf
//	Sinh(NaN) = NaN
public func ndSinh(d: Number): Number {
    if (d.real == 0.0) {
        return Number(
            d.real,
            d.emag,
        )
    }
    if (d.real.isInf()) {
        return Number(
            d.real,
            Float64.Inf,
        )
    }
    let fn = sinh(d.real)
    let deriv = cosh(d.real)
    return Number(
        fn,
        deriv * d.emag,
    )
}

// Cosh returns the hyperbolic cosine of d.
//
// Special cases are:
//
//	Cosh(±0) = 1
//	Cosh(±Inf) = +Inf
//	Cosh(NaN) = NaN
public func ndCosh(d: Number): Number {
    if (d.real.isInf()) {
        return Number(
            Float64.Inf,
            d.real,
        )
    }
    let fn = cosh(d.real)
    let deriv = sinh(d.real)
    return Number(
        fn,
        deriv * d.emag,
    )
}

// Tanh returns the hyperbolic tangent of d.
//
// Special cases are:
//
//	Tanh(±0) = (±0+Nϵ)
//	Tanh(±Inf) = (±1+0ϵ)
//	Tanh(NaN) = NaN
public func ndTanh(d: Number): Number {
    if (d.real == 0.0) {
        return Number(
            d.real,
            d.emag
        )
    } else if (d.real == Float64.Inf) {
        return Number(
            1.0,
            0.0
        )
    } else if (d.real == -Float64.Inf) {
        return Number(
            -1.0,
            0.0
        )
    }
    let fn = tanh(d.real)
    let deriv = 1.0 - fn * fn
    return Number(
        fn,
        deriv * d.emag,
    )
}

// Asinh returns the inverse hyperbolic sine of d.
//
// Special cases are:
//
//	Asinh(±0) = (±0+Nϵ)
//	Asinh(±Inf) = ±Inf
//	Asinh(NaN) = NaN
public func ndAsinh(d: Number): Number {
    if (d.real == 0.0) {
        return Number(
            d.real,
            d.emag
        )
    }
    let fn = asinh(d.real)
    let deriv = 1.0 / sqrt(d.real * d.real + 1.0)
    return Number(
        fn,
        deriv * d.emag,
    )
}

// Acosh returns the inverse hyperbolic cosine of d.
//
// Special cases are:
//
//	Acosh(+Inf) = +Inf
//	Acosh(1) = (0+Infϵ)
//	Acosh(x) = NaN if x < 1
//	Acosh(NaN) = NaN
public func ndAcosh(d: Number): Number {
    if (d.real <= 1.0) {
        if (d.real == 1.0) {
            return Number(
                0.0,
                Float64.Inf,
            )
        }
        return Number(
            Float64.NaN,
            Float64.NaN,
        )
    }
    let fn = acosh(d.real)
    let deriv = 1.0 / sqrt(d.real * d.real - 1.0)
    return Number(
        fn,
        deriv * d.emag,
    )
}

/*
因为仓颉的math库函数atanh不能对(-1,1)以外的数据进行处理，会直接提示异常退出程序影响使用
所以这里改为调用nanAtanh函数，可以实现处理超范围的数据，能对(-1,1)以外的数据返回Float64.NaN
 */
// Atanh returns the inverse hyperbolic tangent of d.
//
// Special cases are:
//
//	Atanh(1) = +Inf
//	Atanh(±0) = (±0+Nϵ)
//	Atanh(-1) = -Inf
//	Atanh(x) = NaN if x < -1 or x > 1
//	Atanh(NaN) = NaN
public func ndAtanh(d: Number): Number {
    if (d.real == 0.0) {
        return Number(
            d.real,
            d.emag,
        )
    }
    if (abs(d.real) == 1.0) {
        if (d.real > 0.0) {
            return Number(
                Float64.Inf,
                Float64.NaN,
            )
        } else {
            return Number(
                -Float64.Inf,
                Float64.NaN,
            )
        }
    }
    let fn = nanAtanh(d.real)
    let deriv = 1.0 / (1.0 - d.real * d.real)
    return Number(
        fn,
        deriv * d.emag,
    )
}
