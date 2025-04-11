package cjnum.cmplxs.cscalar

import std.math.*
import cjnum.floats.scalar
import cjnum.complex.*

//TODO: Cangjie dones't have function to calculating the modulus of a complex128 number
func complexAbs(a: Complex128): Float64 {
    return sqrt(pow(a.imag, 2) + pow(a.real, 2))
}

public func equalWithinAbs(a: Complex128, b: Complex128, tol: Float64): Bool {
    return a == b || complexAbs(Complex128((a.real - b.real), (a.imag- b.imag))) <= tol
}

const minNormalFloat64 = 0x1p-1022

public func equalWithinRel(a: Complex128, b: Complex128, tol: Float64): Bool {
    if (a == b) {
        return true
    }

    let delta = complexAbs(Complex128((a.real - b.real), (a.imag- b.imag))) 

    if (delta <= minNormalFloat64) {
		return delta <= tol * minNormalFloat64
	}
    return delta/max(complexAbs(a), complexAbs(b)) <= tol
}

public func equalWithinAbsOrRel(a: Complex128, b: Complex128, absTol: Float64, reTol: Float64): Bool {
    return equalWithinAbs(a, b, absTol) || equalWithinRel(a, b, reTol)
}

public func parseWithNA(s: String, missing: String): (Complex128, Float64) {
    if (s == missing) {
        return (Complex128(0.0), 0.0)
    }
    let gerParse = cParse(s)
    var weight = 1.0
    return (gerParse, weight)
}

public func round(x: Complex128, prec: Int64): Complex128 {
    if (x == Complex128(0.0)) {
        return Complex128(0.0) 
    }
    return Complex128(scalar.round(x.real, prec), scalar.round(x.imag, prec))
} 

public func roundEven(x: Complex128, prec: Int64): Complex128 {
    if (x == Complex128(0.0)) {
        return Complex128(0.0)
    }
    return Complex128(scalar.roundEven(x.real, prec), scalar.roundEven(x.imag, prec))
}

public func same(a:Complex128, b:Complex128): Bool {
    return (a == b) || ((a.real.isNaN() && a.imag.isNaN()) &&
            (b.real.isNaN() && b.imag.isNaN()))
}