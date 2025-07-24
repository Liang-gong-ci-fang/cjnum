package cjnum.util

import std.math.*

/*
因为函数atan与仓颉库函数atan命名一致，所以这里使用uAtan表示,u表示包名util
 */
// Atan returns the arctangent, in radians, of x.
//
// Special cases are:
//
//	Atan(±0) = ±0
//	Atan(±Inf) = ±Pi/2
public func uAtan(x: Float64): Float64 {
    // TODO if haveArchAtan2 {
    return _atan(x)
}

private func _atan(x: Float64): Float64 {
    if (x == 0.0) {
        return x
    }
    if (x > 0.0) {
        return satan(x)
    }
    return -satan(-x)
}

// Atan2 returns the arc tangent of y/x, using
// the signs of the two to determine the quadrant
// of the return value.
//
// Special cases are (in order):
//
//	Atan2(y, NaN) = NaN
//	Atan2(NaN, x) = NaN
//	Atan2(+0, x>=0) = +0
//	Atan2(-0, x>=0) = -0
//	Atan2(+0, x<=-0) = +Float64.getPI()
//	Atan2(-0, x<=-0) = -Float64.getPI()
//	Atan2(y>0, 0) = +Float64.getPI()/2
//	Atan2(y<0, 0) = -Float64.getPI()/2
//	Atan2(+Inf, +Inf) = +Float64.getPI()/4
//	Atan2(-Inf, +Inf) = -Float64.getPI()/4
//	Atan2(+Inf, -Inf) = 3Pi/4
//	Atan2(-Inf, -Inf) = -3Pi/4
//	Atan2(y, +Inf) = 0
//	Atan2(y>0, -Inf) = +Float64.getPI()
//	Atan2(y<0, -Inf) = -Float64.getPI()
//	Atan2(+Inf, x) = +Float64.getPI()/2
//	Atan2(-Inf, x) = -Float64.getPI()/2

// TODO 100-dev
// public func atan2(y: Float64, x: Float64): Float64 {
//     return _atan2(y, x)
// }

// private func _atan2(y: Float64, x: Float64): Float64 {
//     // special cases
//     if (y.isNaN() || x.isNaN()) {
//         return Float64.NaN
//     } else if (y == 0.0) {
//         if (x >= 0.0 && (!signbit(x))) {
//             return copysign(0.0, y)
//         }
//         return copysign(Float64.getPI(), y)
//     } else if (x == 0.0) {
//         return copysign(Float64.getPI() / 2.0, y)
//     } else if (x.isInf()) {
//         if (x > 0.0) {
//             if (y.isInf()) {
//                 return copysign(Float64.getPI() / 4.0, y)
//             } else {
//                 return copysign(0.0, y)
//             }
//         }
//         if (y.isInf()) {
//             return copysign(Float64.getPI() * 3.0 / 4.0, y)
//         } else {
//             return copysign(Float64.getPI(), y)
//         }
//     } else if (y.isInf()) {
//         return copysign(Float64.getPI() / 2.0, y)
//     }

//     // Call atan and determine the quadrant.
//     // TODO
//     let q = uAtan(y / x)
//     if (x < 0.0) {
//         if (q <= 0.0) {
//             return q + Float64.getPI()
//         }
//         return q - Float64.getPI()
//     }
//     return q
// }

public func copysign(f: Float64, sign: Float64): Float64 {
    let signBit: UInt64 = 1 << 63
    return Float64.fromBits(f.toBits() & (!signBit) | sign.toBits() & signBit)
}

public func nanAcos(x: Float64): Float64 {
    if (x <= 1.0 && x >= -1.0) {
        return acos(x)
    }
    return Float64.NaN
}

public func nanAcosh(x: Float64): Float64 {
    if (x >= 1.0) {
        return acosh(x)
    }
    return Float64.NaN
}

/**
因为仓颉math库函数asin处理在[-1,1]范围外的数的方式是直接抛出异常终止程序，
不利于程序运行，所以这里新建的nanAsin函数能对输入值进行判断，
如果输入值超出[-1,1]范围就返回NaN，否则使用math库函数asin正常处理范围内的数
 */
public func nanAsin(x: Float64): Float64 {
    if (x == 0.0) {
        return x
    }
    if (x <= 1.0 && x >= -1.0) {
        return asin(x)
    }
    return Float64.NaN
}

public func nanAtanh(x: Float64): Float64 {
    if (x == -1.0) {
        return -Float64.Inf
    } else if (x == 1.0) {
        return Float64.Inf
    } else if ((x > 1.0) || (-1.0 > x)) {
        return Float64.NaN
    }
    return atanh(x)
}

public func nanTanh(x: Float64): Float64 {
    if (x == Float64.NaN) {
        return Float64.NaN
    }
    return tanh(x)
}

/*
为了处理仓颉库函数不能处理自变量为Float64.NaN类型，所以这里新建了nanSinh函数能够处理
自变量为Float64.NaN类型的情况
*/
public func nanSinh(x: Float64): Float64 {
    if (x == Float64.NaN) {
        return Float64.NaN
    }
    return sinh(x)
}

/*
为了处理仓颉库函数不能处理自变量为Float64.NaN类型，所以这里新建了nanSinh函数能够处理
自变量为Float64.NaN类型的情况
*/
public func nanCosh(x: Float64): Float64 {
    if (x == Float64.NaN) {
        return Float64.NaN
    }
    return cosh(x)
}

/**
因为仓颉math库函数sqrt处理小于0的数的方式是直接抛出异常终止程序，
不利于程序运行，所以这里新建的nanAsin函数能对输入值进行判断，
如果输入值小于0就返回NaN，否则使用math库函数asin正常处理范围内的数
 */
public func nanSqrt(x: Float64): Float64 {
    if (x < 0.0) {
        return Float64.NaN
    }
    return sqrt(x)
}

// satan reduces its argument (known to be positive)
// to the range [0, 0.66] and calls xatan.
private func satan(x: Float64): Float64 {
    let morebits = 6.123233995736765886130e-17 // pi/2 = PIO2 + morebits
    let tan3pio8 = 2.41421356237309504880 // tan(3*pi/8)
    if (x <= 0.66) {
        return xatan(x)
    }
    if (x > tan3pio8) {
        return Float64.getPI() / 2.0 - xatan(1.0 / x) + morebits
    }
    return Float64.getPI() / 4.0 + xatan((x - 1.0) / (x + 1.0)) + 0.5 * morebits
}

public func signbit(x: Float64): Bool {
    return (x.toBits() & (1 << 63)) != 0
}

let reduceThreshold = 1 << 29

// Coefficients _sin[] and _cos[] are found in pkg/math/sin.go.

// Sincos returns Sin(x), Cos(x).
//
// Special cases are:
//
//	Sincos(±0) = ±0, 1
//	Sincos(±Inf) = NaN, NaN
//	Sincos(NaN) = NaN, NaN
public func sincos(x_: Float64): (sinX: Float64, cosX: Float64) {
    /*
    因为仓颉语言为x_不可修改类型，所以这里创建可修改变量x
     */
    var x = x_
    let PI4A = 7.85398125648498535156e-1 // 0x3fe921fb40000000, Float64.getPI()/4 split into three parts
    let PI4B = 3.77489470793079817668e-8 // 0x3e64442d00000000,
    let PI4C = 2.69515142907905952645e-15 // 0x3ce8469898cc5170,

    // special cases
    if (x == 0.0) {
        return (x, 1.0) // return ±0.0, 1.0
    } else if (x.isNaN() || x.isInf()) {
        return (Float64.NaN, Float64.NaN)
    }

    // make argument positive
    var sinSign = false
    var cosSign = false
    if (x < 0.0) {
        x = -x
        sinSign = true
    }

    var j: UInt64 = 0
    var y: Float64 = 0.0
    var z: Float64 = 0.0
    if (x >= Float64(reduceThreshold)) {
        (j, z) = trigReduce(x)
    } else {
        j = UInt64(x * (4.0 / Float64.getPI())) // integer part of x/(Float64.getPI()/4), as integer for tests on the phase angle
        y = Float64(j) // integer part of x/(Float64.getPI()/4), as float

        if ((j & 1) == 1) { // map zeros to origin
            j++
            y += 1.0
        }
        j &= 7 // octant modulo 2Pi radians (360 degrees)
        z = ((x - y * PI4A) - y * PI4B) - y * PI4C // Extended precision modular arithmetic
    }
    if (j > 3) { // reflect in x axis
        j -= 4
        sinSign = !sinSign
        cosSign = !cosSign
    }
    if (j > 1) {
        cosSign = !cosSign
    }

    let zz = z * z
    var cos = 1.0 - 0.5 * zz + zz * zz * ((((((_cos[0] * zz) + _cos[1]) * zz + _cos[2]) * zz + _cos[3]) * zz + _cos[4]) *
        zz + _cos[5])
    var sin = z + z * zz * ((((((_sin[0] * zz) + _sin[1]) * zz + _sin[2]) * zz + _sin[3]) * zz + _sin[4]) * zz + _sin[5])
    if (j == 1 || j == 2) {
        var temp = sin
        sin = cos
        cos = temp
    }
    if (cosSign) {
        cos = -cos
    }
    if (sinSign) {
        sin = -sin
    }
    return (sin, cos)
}

// calculate sinh and cosh.
public func sinhcosh(x: Float64): (Float64, Float64) {
    if (abs(x) <= 0.5) {
        return (sinh(x), cosh(x))
    }
    var e = exp(x)
    let ei = 0.5 / e
    e *= 0.5
    return (e - ei, e + ei)
}

// mPi4 is the binary digits of 4/pi as a uint64 array,
// that is, 4/pi = Sum mPi4[i]*2^(-64*i)
// 19 64-bit digits and the leading one bit give 1217 bits
// of precision to handle the largest possible float64 exponent.
let mPi4: Array<UInt64> = [
    0x0000000000000001,
    0x45f306dc9c882a53,
    0xf84eafa3ea69bb81,
    0xb6c52b3278872083,
    0xfca2c757bd778ac3,
    0x6e48dc74849ba5c0,
    0x0c925dd413a32439,
    0xfc3bd63962534e7d,
    0xd1046bea5d768909,
    0xd338e04d68befc82,
    0x7323ac7306a673e9,
    0x3908bf177bf25076,
    0x3ff12fffbc0b301f,
    0xde5e2316b414da3e,
    0xda6cfd9e4f96136e,
    0x9e8c7ecd3cbfd45a,
    0xea4f758fd7cbe2f6,
    0x7a0e73ef14a525d4,
    0xd7f6bf623f1aba10,
    0xac06608df8f6d757
]

// trigReduce implements Payne-Hanek range reduction by Float64.getPI()/4
// for x > 0. It returns the integer part mod 8 (j) and
// the fractional part (z) of x / (Float64.getPI()/4).
// The implementation is based on:
// "ARGUMENT REDUCTION FOR HUGE ARGUMENTS: Good to the Last Bit"
// K. C. Ng et al, March 24, 1992
// The simulated multi-precision calculation of x*B uses 64-bit integer arithmetic.
public func trigReduce(x: Float64): (j: UInt64, z: Float64) {
    let PI4 = Float64.getPI() / 4.0
    if (x < PI4) {
        return (0, x)
    }
    // Extract out the integer and exponent such that,
    // x = ix * 2 ** exp.
    var ix = x.toBits()
    let exp = Int64(ix >> shift & mask - bias - shift)
    ix = ix & UInt64(!(mask << shift))
    ix |= 1 << shift
    // Use the exponent to extract the 3 appropriate uint64 digits from mPi4,
    // B ~ (z0, z1, z2), such that the product leading digit has the exponent -61.
    // Note, exp >= -53 since x >= PI4 and exp < 971 for maximum float64.
    let digit = Int64(exp + 61) / 64
    let bitshift = UInt64(exp + 61) % 64
    let z0 = (mPi4[digit] << bitshift) | (mPi4[digit + 1] >> (64 - bitshift))
    let z1 = (mPi4[digit + 1] << bitshift) | (mPi4[digit + 2] >> (64 - bitshift))
    let z2 = (mPi4[digit + 2] << bitshift) | (mPi4[digit + 3] >> (64 - bitshift))
    // Multiply mantissa by the digits and extract the upper two digits (hi, lo).
    let (z2hi, _) = mul64(z2, ix)
    let (z1hi, z1lo) = mul64(z1, ix)
    let z0lo = z0 * ix
    let (lo, c) = add64(z1lo, z2hi, 0)
    var (hi, _) = add64(z0lo, z1hi, c)
    // The top 3 bits are j.
    var j = hi >> 61
    // Extract the fraction and find its magnitude.
    hi = hi << 3 | lo >> 61
    let lz = UInt64(leadingZeros64(hi))
    let e = UInt64(bias - (lz + 1))
    // Clear implicit mantissa bit and shift into place.
    hi = (hi << (lz + 1)) | (lo >> (64 - (lz + 1)))
    hi >>= 64 - shift
    // Include the exponent and convert to a float.
    hi |= e << shift
    var z = Float64.fromBits(hi)
    // Map zeros to origin.
    if ((j & 1) == 1) {
        j++
        j &= 7
        z -= 1.0
    }
    // Multiply the fractional part by pi/4.
    return (j, z * PI4)
}

// xatan evaluates a series valid in the range [0, 0.66].
private func xatan(x: Float64): Float64 {
    let P0 = -8.750608600031904122785e-1
    let P1 = -1.615753718733365076637e1
    let P2 = -7.500855792314704667340e1
    let P3 = -1.228866684490136173410e2
    let P4 = -6.485021904942025371773e1
    let Q0 = 2.485846490142306297962e1
    let Q1 = 1.650270098316988542046e2
    let Q2 = 4.328810604912902668951e2
    let Q3 = 4.853903996359136964868e2
    let Q4 = 1.945506571482613964425e2

    var z = x * x
    z = z * ((((P0 * z + P1) * z + P2) * z + P3) * z + P4) / (((((z + Q0) * z + Q1) * z + Q2) * z + Q3) * z + Q4)
    z = x * z + x
    return z
}
