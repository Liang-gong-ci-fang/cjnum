package cjnum.util

import std.math.*

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
