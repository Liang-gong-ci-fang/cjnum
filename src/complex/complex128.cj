package cjnum.complex

import cjnum.util.*
import std.math.*

public class Complex128 <: ToString {
    public static let NaN: Complex128 = Complex128(Float64.NaN, Float64.NaN)
    public static let Inf: Complex128 = Complex128(Float64.Inf, Float64.Inf)

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
            (this.real * right),
            (this.imag * right)
        )
    }

    public operator func /(right: Float64) {
        if (right == 0.0) {
            return Complex128(
                Float64.NaN,
                Float64.NaN
            )
        }
        return Complex128(
            (this.real / right),
            (this.imag / right)
        )
    }

    public func toString(): String {
        if (this.isNaN()) {
            return Float64.NaN.toString()
        }
        if (this.isInf()) {
            return Float64.Inf.toString()
        }

        let r = this.real
        let i = this.imag

        if (i == 0.0) {
            return r.toString()
        }

        if (r == 0.0) {
            return if (i == 1.0) {
                "i"
            } else if (i == -1.0) {
                "-i"
            } else {
                i.toString() + "i"
            }
        }

        let sign = if (i > 0.0) {
            "+"
        } else {
            "-"
        }
        let absI = abs(i)

        if (absI == 1.0) {
            return r.toString() + sign + "i"
        }

        return r.toString() + sign + absI.toString() + "i"
    }

    // TODO: add more operators and functions
}

// Abs returns the absolute value (also called the modulus) of x.
public func c128Abs(x: Complex128): Float64 {
    return hypot(x.real, x.imag)
}

// Asin returns the inverse sine of x.
public func c128Asin(x: Complex128): Complex128 {
    let re = x.real
    let im = x.imag
    if (im == 0.0 && abs(re) <= 1.0) {
        return Complex128(asin(re), im)
    } else if (re == 0.0 && abs(im) <= 1.0) {
        return Complex128(re, asinh(im))
    } else if (im.isNaN()) {
        if (re == 0.0) {
            return Complex128(re, Float64.NaN)
        } else if (re.isInf()) {
            return Complex128(Float64.NaN, re)
        } else {
            return Complex128.NaN
        }
    } else if (im.isInf()) {
        if (re.isNaN()) {
            return x
        } else if (re.isInf()) {
            return Complex128(copysign(Float64.PI / 4.0, re), im)
        } else {
            return Complex128(copysign(0.0, re), im)
        }
    } else if (re.isInf()) {
        return Complex128(copysign(Float64.PI / 2.0, re), copysign(re, im))
    }
    let ct = Complex128(-x.imag, x.real) // i * x
    let xx = x * x
    let x1 = Complex128(1.0 - xx.real, -xx.imag) // 1 - x*x
    let x2 = c128Sqrt(x1) // x2 = sqrt(1 - x*x)
    let w = c128Log(ct + x2)
    return Complex128(w.imag, -w.real) // -i * w
}

// Asinh returns the inverse hyperbolic sine of x.
public func c128Asinh(x: Complex128): Complex128 {
    let re = x.real
    let im = x.imag
    if (im == 0.0 && abs(re) <= 1.0) {
        return Complex128(asinh(re), im)
    } else if (re == 0.0 && abs(im) <= 1.0) {
        return Complex128(re, asin(im))
    } else if (re.isInf()) {
        if (im.isInf()) {
            return Complex128(re, copysign(Float64.PI / 4.0, im))
        } else if (im.isNaN()) {
            return x
        } else {
            return Complex128(re, copysign(0.0, im))
        }
    } else if (re.isNaN()) {
        if (im == 0.0) {
            return x
        } else if (im.isInf()) {
            return Complex128(im, re)
        } else {
            return Complex128.NaN
        }
    } else if (im.isInf()) {
        return Complex128(copysign(im, re), copysign(Float64.PI / 2.0, im))
    }
    let xx = x * x
    let x1 = Complex128(1.0 + xx.real, xx.imag) // 1 + x*x
    return c128Log(x + c128Sqrt(x1)) // log(x + sqrt(1 + x*x))
}

// Complex circular arc cosine
//
// DESCRIPTION:
//
// w = arccos z  =  PI/2 - arcsin z.
//
// ACCURACY:
//
//                      Relative error:
// arithmetic   domain     # trials      peak         rms
//    DEC       -10,+10      5200      1.6e-15      2.8e-16
//    IEEE      -10,+10     30000      1.8e-14      2.2e-15

// Acos returns the inverse cosine of x.
public func c128Acos(x: Complex128): Complex128 {
    let w = c128Asin(x)
    return Complex128(Float64.PI / 2.0 - w.real, -w.imag)
}

// Acosh returns the inverse hyperbolic cosine of x.
public func c128Acosh(x: Complex128): Complex128 {
    if (x == Complex128(0.0, 0.0)) {
        return Complex128(0.0, copysign(Float64.PI / 2.0, x.imag))
    }
    let w = c128Acos(x)
    if (w.imag <= 0.0) {
        return Complex128(-w.imag, w.real) // i * w
    }
    return Complex128(w.imag, -w.real) // -i * w
}

// Complex circular arc tangent
//
// DESCRIPTION:
//
// If
//     z = x + iy,
//
// then
//          1       (    2x     )
// Re w  =  - arctan(-----------)  +  k PI
//          2       (     2    2)
//                  (1 - x  - y )
//
//               ( 2         2)
//          1    (x  +  (y+1) )
// Im w  =  - log(------------)
//          4    ( 2         2)
//               (x  +  (y-1) )
//
// Where k is an arbitrary integer.
//
// catan(z) = -i catanh(iz).
//
// ACCURACY:
//
//                      Relative error:
// arithmetic   domain     # trials      peak         rms
//    DEC       -10,+10      5900       1.3e-16     7.8e-18
//    IEEE      -10,+10     30000       2.3e-15     8.5e-17
// The check catan( ctan(z) )  =  z, with |x| and |y| < PI/2,
// had peak relative error 1.5e-16, rms relative error
// 2.9e-17.  See also clog().

// Atan returns the inverse tangent of x.
public func c128Atan(x: Complex128): Complex128 {
    let re = x.real
    let im = x.imag
    if (im == 0.0) {
        return Complex128(atan(re), im)
    } else if (re == 0.0 && abs(im) <= 1.0) {
        return Complex128(re, atanh(im))
    } else if (im.isInf() || re.isInf()) {
        if (re.isNaN()) {
            return Complex128(Float64.NaN, copysign(0.0, im))
        }
        return Complex128(copysign(Float64.PI / 2.0, re), copysign(0.0, im))
    } else if (re.isNaN() || im.isInf()) {
        return Complex128.NaN
    }
    let x2 = x.real * x.real
    let a = 1.0 - x2 - x.imag * x.imag
    if (a == 0.0) {
        return Complex128.NaN
    }
    var t = 0.5 * atan2(2.0 * x.real, a)
    let w = reducePi(t)

    t = x.imag - 1.0
    let b = x2 + t * t
    if (b == 0.0) {
        return Complex128.NaN
    }
    t = x.imag + 1.0
    let c = (x2 + t * t) / b
    return Complex128(w, 0.25 * log(c))
}

// Atanh returns the inverse hyperbolic tangent of x.
public func c128Atanh(x: Complex128): Complex128 {
    var z = Complex128(-x.imag, x.real) // z = i * x
    z = c128Atan(z)
    return Complex128(z.imag, -z.real) // z = -i * z
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
    let re = x.real
    let im = x.imag
    if (re.isInf()) {
        if (re > 0.0 && im == 0.0) {
            return x
        } else if (im.isInf() || im.isNaN()) {
            if (re < 0.0) {
                return Complex128(0.0, copysign(0.0, im))
            } else {
                return Complex128(Float64.Inf, Float64.NaN)
            }
        }
    } else if (re.isNaN()) {
        if (im == 0.0) {
            return Complex128(Float64.NaN, im)
        }
    }
    let r = exp(x.real)
    let (s, c) = sincos(x.imag)
    return Complex128(r * c, r * s)
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

// Sin returns the sine of x.
public func c128Sin(x: Complex128): Complex128 {
    let re = x.real
    let im = x.imag
    if (im == 0.0 && (re.isInf() || re.isNaN())) {
        return Complex128(Float64.NaN, im)
    } else if (im.isInf()) {
        if (re == 0.0) {
            return x
        } else if (re.isInf() || re.isNaN()) {
            return Complex128(Float64.NaN, im)
        }
    } else if (re == 0.0 && im.isNaN()) {
        return x
    }
    let (s, c) = sincos(x.real)
    let (sh, ch) = sinhcosh(x.imag)
    return Complex128(s * ch, c * sh)
}

// Complex hyperbolic sine
//
// DESCRIPTION:
//
// csinh z = (cexp(z) - cexp(-z))/2
//         = sinh x * cos y  +  i cosh x * sin y .
//
// ACCURACY:
//
//                      Relative error:
// arithmetic   domain     # trials      peak         rms
//    IEEE      -10,+10     30000       3.1e-16     8.2e-17

// Sinh returns the hyperbolic sine of x.
public func c128Sinh(x: Complex128): Complex128 {
    let re = x.real
    let im = x.imag
    if (re == 0.0 && (im.isInf() || im.isNaN())) {
        return Complex128(re, Float64.NaN)
    } else if (re.isInf()) {
        if (im == 0.0) {
            return Complex128(re, im)
        } else if (im.isInf() || im.isNaN()) {
            return Complex128(re, Float64.NaN)
        }
    } else if (im == 0.0 && re.isNaN()) {
        return Complex128(Float64.NaN, im)
    }
    let (s, c) = sincos(x.imag)
    let (sh, ch) = sinhcosh(x.real)
    return Complex128(c * sh, s * ch)
}

// Complex circular cosine
//
// DESCRIPTION:
//
// If
//     z = x + iy,
//
// then
//
//     w = cos x  cosh y  -  i sin x sinh y.
//
// ACCURACY:
//
//                      Relative error:
// arithmetic   domain     # trials      peak         rms
//    DEC       -10,+10      8400       4.5e-17     1.3e-17
//    IEEE      -10,+10     30000       3.8e-16     1.0e-16

// Cos returns the cosine of x.
public func c128Cos(x: Complex128): Complex128 {
    let re = x.real
    let im = x.imag
    if (im == 0.0 && (re.isInf() || re.isNaN())) {
        return Complex128(Float64.NaN, -im * copysign(0.0, re))
    } else if (im.isInf()) {
        if (re == 0.0) {
            return Complex128(Float64.Inf, -re * copysign(0.0, im))
        } else if (re.isInf() || re.isNaN()) {
            return Complex128(Float64.Inf, Float64.NaN)
        }
    } else if (re == 0.0 && im.isNaN()) {
        return Complex128(Float64.NaN, 0.0)
    }
    let (s, c) = sincos(x.real)
    let (sh, ch) = sinhcosh(x.imag)
    return Complex128(c * ch, -s * sh)
}

// Complex hyperbolic cosine
//
// DESCRIPTION:
//
// ccosh(z) = cosh x  cos y + i sinh x sin y .
//
// ACCURACY:
//
//                      Relative error:
// arithmetic   domain     # trials      peak         rms
//    IEEE      -10,+10     30000       2.9e-16     8.1e-17

// Cosh returns the hyperbolic cosine of x.
public func c128Cosh(x: Complex128): Complex128 {
    let re = x.real
    let im = x.imag
    if (re == 0.0 && (im.isInf() || im.isNaN())) {
        return Complex128(Float64.NaN, re * copysign(0.0, im))
    } else if (re.isInf()) {
        if (im == 0.0) {
            return Complex128(Float64.Inf, im * copysign(0.0, re))
        } else if (im.isInf() || im.isNaN()) {
            return Complex128(Float64.Inf, Float64.NaN)
        }
    } else if (im == 0.0 && re.isNaN()) {
        return Complex128(Float64.NaN, im)
    }
    let (s, c) = sincos(x.imag)
    let (sh, ch) = sinhcosh(x.real)
    return Complex128(c * ch, s * sh)
}

// Complex circular tangent
//
// DESCRIPTION:
//
// If
//     z = x + iy,
//
// then
//
//           sin 2x  +  i sinh 2y
//     w  =  --------------------.
//            cos 2x  +  cosh 2y
//
// On the real axis the denominator is zero at odd multiples
// of PI/2. The denominator is evaluated by its Taylor
// series near these points.
//
// ctan(z) = -i ctanh(iz).
//
// ACCURACY:
//
//                      Relative error:
// arithmetic   domain     # trials      peak         rms
//    DEC       -10,+10      5200       7.1e-17     1.6e-17
//    IEEE      -10,+10     30000       7.2e-16     1.2e-16
// Also tested by ctan * ccot = 1 and catan(ctan(z))  =  z.

// Tan returns the tangent of x.
public func c128Tan(x: Complex128): Complex128 {
    let re = x.real
    let im = x.imag
    if (im.isInf()) {
        if (re.isInf() || re.isNaN()) {
            return Complex128(copysign(0.0, re), copysign(1.0, im))
        }
        return Complex128(copysign(0.0, sin(2.0 * re)), copysign(1.0, im))
    } else if (re == 0.0 && im.isNaN()) {
        return x
    }
    var d = cos(2.0 * x.real) + cosh(2.0 * x.imag)
    if (abs(d) < 0.25) {
        d = tanSeries(x)
    }
    if (d == 0.0) {
        return Complex128.Inf
    }
    return Complex128(sin(2.0 * x.real) / d, sinh(2.0 * x.imag) / d)
}

// Complex hyperbolic tangent
//
// DESCRIPTION:
//
// tanh z = (sinh 2x  +  i sin 2y) / (cosh 2x + cos 2y) .
//
// ACCURACY:
//
//                      Relative error:
// arithmetic   domain     # trials      peak         rms
//    IEEE      -10,+10     30000       1.7e-14     2.4e-16

// Tanh returns the hyperbolic tangent of x.
public func c128Tanh(x: Complex128): Complex128 {
    let re = x.real
    let im = x.imag
    if (re.isInf()) {
        if (im.isInf() || im.isNaN()) {
            return Complex128(copysign(1.0, re), copysign(0.0, im))
        }
        return Complex128(copysign(1.0, re), copysign(0.0, sin(2.0 * im)))
    } else if (im == 0.0 && re.isNaN()) {
        return x
    }
    var d = cos(2.0 * x.real) + cosh(2.0 * x.imag)
    if (d == 0.0) {
        return Complex128.Inf
    }
    return Complex128(sinh(2.0 * x.real) / d, sin(2.0 * x.imag) / d)
}

// reducePi reduces the input argument x to the range (-Pi/2, Pi/2].
// x must be greater than or equal to 0. For small arguments it
// uses Cody-Waite reduction in 3 float64 parts based on:
// "Elementary Function Evaluation:  Algorithms and Implementation"
// Jean-Michel Muller, 1997.
// For very large arguments it uses Payne-Hanek range reduction based on:
// "ARGUMENT REDUCTION FOR HUGE ARGUMENTS: Good to the Last Bit"
// K. C. Ng et al, March 24, 1992.
private func reducePi(x_: Float64): Float64 {
    var x = x_
    // reduceThreshold is the maximum value of x where the reduction using
    // Cody-Waite reduction still gives accurate results. This threshold
    // is set by t*PIn being representable as a float64 without error
    // where t is given by t = floor(x * (1 / Pi)) and PIn are the leading partial
    // terms of Pi. Since the leading terms, PI1 and PI2 below, have 30 and 32
    // trailing zero bits respectively, t should have less than 30 significant bits.
    //	t < 1<<30  -> floor(x*(1/Pi)+0.5) < 1<<30 -> x < (1<<30-1) * Pi - 0.5
    // So, conservatively we can take x < 1<<30.
    let reduceThreshold: Float64 = Float64(1 << 30)
    if (abs(x) < reduceThreshold) {
        // Use Cody-Waite reduction in three parts.
        // PI1, PI2 and PI3 comprise an extended precision value of PI
        // such that PI ~= PI1 + PI2 + PI3. The parts are chosen so
        // that PI1 and PI2 have an approximately equal number of trailing
        // zero bits. This ensures that t*PI1 and t*PI2 are exact for
        // large integer values of t. The full precision PI3 ensures the
        // approximation of PI is accurate to 102 bits to handle cancellation
        // during subtraction.
        let PI1 = 3.141592502593994 // 0x400921fb40000000
        let PI2 = 1.5099578831723193e-07 // 0x3e84442d00000000
        let PI3 = 1.0780605716316238e-14 // 0x3d08469898cc5170

        var t = x / Float64.PI
        t += 0.5
        t = Float64(Int64(t)) // int64(t) = the multiple
        return ((x - t * PI1) - t * PI2) - t * PI3
    }
    // Must apply Payne-Hanek range reduction
    let mask: UInt64 = 0x7FF
    let shift: UInt64 = 64 - 11 - 1
    let bias: UInt64 = 1023
    let fracMask: UInt64 = 1 << shift - 1

    // Extract out the integer and exponent such that,
    // x = ix * 2 ** exp.
    var ix = x.toBits()
    let exp: Int64 = Int64(ix >> shift & mask - bias - shift)
    ix &= fracMask
    ix |= 1 << shift

    // mPi is the binary digits of 1/Pi as a uint64 array,
    // that is, 1/Pi = Sum mPi[i]*2^(-64*i).
    // 19 64-bit digits give 1216 bits of precision
    // to handle the largest possible float64 exponent.
    let mPi: Array<UInt64> = [
        0x0000000000000000,
        0x517cc1b727220a94,
        0xfe13abe8fa9a6ee0,
        0x6db14acc9e21c820,
        0xff28b1d5ef5de2b0,
        0xdb92371d2126e970,
        0x0324977504e8c90e,
        0x7f0ef58e5894d39f,
        0x74411afa975da242,
        0x74ce38135a2fbf20,
        0x9cc8eb1cc1a99cfa,
        0x4e422fc5defc941d,
        0x8ffc4bffef02cc07,
        0xf79788c5ad05368f,
        0xb69b3f6793e584db,
        0xa7a31fb34f2ff516,
        0xba93dd63f5f2f8bd,
        0x9e839cfbc5294975,
        0x35fdafd88fc6ae84,
        0x2b0198237e3db5d5
    ]
    // Use the exponent to extract the 3 appropriate uint64 digits from mPi,
    // B ~ (z0, z1, z2), such that the product leading digit has the exponent -64.
    // Note, exp >= 50 since x >= reduceThreshold and exp < 971 for maximum float64.
    let digit = UInt64(exp + 64) / 64
    let bitshift = UInt64(exp + 64) % 64
    let z0 = (mPi[Int64(digit)] << bitshift) | (mPi[Int64(digit) + 1] >> (64 - bitshift))
    let z1 = (mPi[Int64(digit) + 1] << bitshift) | (mPi[Int64(digit) + 2] >> (64 - bitshift))
    let z2 = (mPi[Int64(digit) + 2] << bitshift) | (mPi[Int64(digit) + 3] >> (64 - bitshift))
    // Multiply mantissa by the digits and extract the upper two digits (hi, lo).
    let (z2hi, _) = mul64(z2, ix)
    let (z1hi, z1lo) = mul64(z1, ix)
    let z0lo = z0 * ix
    let (lo, c) = add64(z1lo, z2hi, 0)
    var (hi, _) = add64(z0lo, z1hi, c)
    // Find the magnitude of the fraction.
    let lz = UInt64(leadingZeros64(hi))
    let e = UInt64(bias - (lz + 1))
    // Clear implicit mantissa bit and shift into place.
    hi = (hi << (lz + 1)) | (lo >> (64 - (lz + 1)))
    hi >>= 64 - shift
    // Include the exponent and convert to a float.
    hi |= e << shift
    x = Float64.fromBits(hi)
    // map to (-Pi/2, Pi/2]
    if (x > 0.5) {
        x -= 1.0
    }
    return Float64.PI * x
}

// Taylor series expansion for cosh(2y) - cos(2x)
private func tanSeries(z: Complex128): Float64 {
    let MACHEP = 1.0 / Float64(1 << 53)
    var x = abs(2.0 * z.real)
    var y = abs(2.0 * z.imag)
    x = reducePi(x)
    x = x * x
    y = y * y
    var x2 = 1.0
    var y2 = 1.0
    var f = 1.0
    var rn = 0.0
    var d = 0.0
    while (true) {
        rn += 1.0
        f *= rn
        rn += 1.0
        f *= rn
        x2 *= x
        y2 *= y
        var t = y2 + x2
        t /= f
        d += t

        rn += 1.0
        f *= rn
        rn += 1.0
        f *= rn
        x2 *= x
        y2 *= y
        t = y2 - x2
        t /= f
        d += t
        if (!(abs(t / d) > MACHEP)) {
            // Caution: Use ! and > instead of <= for correct behavior if t/d is NaN.
            // See issue 17577.
            break
        }
    }
    return d
}

// Complex square root
//
// DESCRIPTION:
//
// If z = x + iy,  r = |z|, then
//
//                       1/2
// Re w  =  [ (r + x)/2 ]   ,
//
//                       1/2
// Im w  =  [ (r - x)/2 ]   .
//
// Cancellation error in r-x or r+x is avoided by using the
// identity  2 Re w Im w  =  y.
//
// Note that -w is also a square root of z. The root chosen
// is always in the right half plane and Im w has the same sign as y.
//
// ACCURACY:
//
//                      Relative error:
// arithmetic   domain     # trials      peak         rms
//    DEC       -10,+10     25000       3.2e-17     9.6e-18
//    IEEE      -10,+10   1,000,000     2.9e-16     6.1e-17

// Sqrt returns the square root of x.
// The result r is chosen so that real(r) ≥ 0 and imag(r) has the same sign as x.imag.
public func c128Sqrt(x: Complex128): Complex128 {
    if (x.imag == 0.0) {
        // Ensure that imag(r) has the same sign as x.imag for x.imag == signed zero.
        if (x.real == 0.0) {
            return Complex128(0.0, x.imag)
        }
        if (x.real < 0.0) {
            return Complex128(0.0, copysign(sqrt(-x.real), x.imag))
        }
        return Complex128(sqrt(x.real), x.imag)
    } else if (x.imag.isInf()) {
        return Complex128(Float64.Inf, x.imag)
    }
    if (x.real == 0.0) {
        if (x.imag < 0.0) {
            let r = sqrt(-0.5 * x.imag)
            return Complex128(r, -r)
        }
        let r = sqrt(0.5 * x.imag)
        return Complex128(r, r)
    }
    var a = x.real
    var b = x.imag
    var scale: Float64 = 0.0
    // Rescale to avoid internal overflow or underflow.
    if (abs(a) > 4.0 || abs(b) > 4.0) {
        a *= 0.25
        b *= 0.25
        scale = 2.0
    } else {
        a *= 1.8014398509481984e16 // 2**54
        b *= 1.8014398509481984e16
        scale = 7.450580596923828125e-9 // 2**-27
    }
    var r = hypot(a, b)
    var t: Float64 = 0.0
    if (a > 0.0) {
        t = sqrt(0.5 * r + 0.5 * a)
        r = scale * abs((0.5 * b) / t)
        t *= scale
    } else {
        r = sqrt(0.5 * r - 0.5 * a)
        t = scale * abs((0.5 * b) / r)
        r *= scale
    }
    if (b < 0.0) {
        return Complex128(t, -r)
    }
    return Complex128(t, r)
}
