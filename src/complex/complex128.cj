package cjnum.complex

import cjnum.util.*
import std.math.*

public class Complex128 {
    public static let NaN: Complex128 = Complex128(Float64.NaN, Float64.NaN)
    public static let Inf: Complex128 = Complex128(Float64.Inf,Float64.Inf)

    public init(real: Float64) {
        this(real, 0.0)
    }

    public Complex128(
        public let real: Float64,
        public let imag: Float64
    ) {}

    public operator func ==(right: Complex128) {
        return this.real == right.real && this.imag == right.imag
    }

    public func isNaN(): Bool {
        return this.real.isNaN() || this.imag.isNaN()
    }
    
    public func isInf(): Bool {
	if (this.real.isInf() || this.imag.isInf()) {
		return true
	}
	    return false
    }

    public func conj(): Complex128 {
        return Complex128(this.real, -this.imag)
    }

    public operator func !=(right: Complex128) {
        return this.real != right.real || this.imag != right.imag
    }

    public operator func -(right: Complex128) {
        return Complex128((this.real - right.real), (this.imag - right.imag))
    }

    public operator func +(right: Complex128) {
        return Complex128((this.real + right.real), (this.imag + right.imag))
    }

    public operator func *(right: Complex128) {
        return Complex128(
            (this.real * right.real - this.imag * right.imag),
            (this.real * right.imag + this.imag * right.real)
        )
    }

    public operator func /(right: Complex128) {
        let denominator = right.real * right.real + right.imag * right.imag
        return Complex128(
            (this.real * right.real + this.imag * right.imag) / denominator,
            (this.imag * right.real - this.real * right.imag) / denominator
        )
    }

    public operator func *(right: Float64) {
        return Complex128(
            (this.real * right) ,
            (this.imag * right)
        )
    }

    public operator func /(right: Float64) {
        if(right==0.0){
            return Complex128(
                Float64.NaN,
                Float64.NaN
            )
        }
        return Complex128(
            (this.real / right) ,
            (this.imag / right)
        )
    }
    // TODO: add more operators and functions
}

// Abs returns the absolute value (also called the modulus) of x.
public func c128Abs(x: Complex128): Float64 { 
    return hypot(x.real, x.imag) 
}

// Complex exponential function
//
// DESCRIPTION:
//
// Returns the complex exponential of the complex argument z.
//
// If
//     z = x + iy,
//     r = exp(x),
// then
//     w = r cos y + i r sin y.
//
// ACCURACY:
//
//                      Relative error:
// arithmetic   domain     # trials      peak         rms
//    DEC       -10,+10      8700       3.7e-17     1.1e-17
//    IEEE      -10,+10     30000       3.0e-16     8.7e-17

// Exp returns e**x, the base-e exponential of x.
public func c128Exp(x: Complex128): Complex128 {
    let re=x.real
    let im=x.imag
    if(re.isInf()){
        if(re>0.0&&im==0.0){
            return x
        }
        else if(im.isInf()||im.isNaN()){
            if(re<0.0){
                return Complex128(0.0,copysign(0.0, im))
            }else {
                return Complex128(Float64.Inf, Float64.NaN)
            }
        }
    }
    else if(re.isNaN()){
        if(im==0.0){
            return Complex128(Float64.NaN,im)
        }
    }
    let r = exp(x.real)
    let (s, c) = sincos(x.imag)
    return Complex128(r*c, r*s)
}

// Log returns the natural logarithm of x.
public func c128Log(x: Complex128): Complex128 {
	return Complex128(log(c128Abs(x)), c128Phase(x))
}

// Phase returns the phase (also called the argument) of x.
// The returned value is in the range [-Pi, Pi].
public func c128Phase(x: Complex128): Float64 { 
    return atan2(x.imag, x.real)
}
