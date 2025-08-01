package cjnum.internal.asm.f64

import std.math.*

public func l1Norm(x: Array<Float64>): Float64 {
    var sum: Float64 = 0.0
    for (v in x) {
        sum += abs(v)
    }
    return sum
}

public func l1NormInc(x: Array<Float64>, n: Int64, incX: Int64): Float64 {
    var sum: Float64 = 0.0
    for (i in 0..n) {
        sum += abs(x[i * incX])
    }
    return sum
}

public func add(dst: Array<Float64>, s: Array<Float64>): Unit {
    for (i in 0..s.size) {
        dst[i] += s[i]
    }
}

public func addConst(alpha: Float64, x: Array<Float64>): Unit {
    for (i in 0..x.size) {
        x[i] += alpha
    }
}

public func cumSum(dst: Array<Float64>, s: Array<Float64>): Array<Float64> {
    if(s.isEmpty()) { 
        return dst
    }
    dst[0] = s[0]
    for (i in 1..s.size) {
        dst[i] = dst[i - 1] + s[i]
    }
    return dst
}

public func cumProd(dst: Array<Float64>, s: Array<Float64>): Array<Float64> {
    if(s.isEmpty()) { 
        return dst
    }
    dst[0] = s[0]
    for (i in 1..s.size) {
        dst[i] = dst[i - 1] * s[i]
    }
    return dst
}

public func div(dst: Array<Float64>, s: Array<Float64>): Unit {
    for (i in 0..s.size) {
        dst[i] = dst[i] / s[i]
    }
}

public func divTo(dst: Array<Float64>, s: Array<Float64>, t: Array<Float64>): Array<Float64> {
    for (i in 0..s.size) {
        dst[i] = s[i] / t[i]
    }
    return dst
}

public func l1Dist(s: Array<Float64>, t: Array<Float64>): Float64 {
    var norm: Float64 = 0.0
    for (i in 0..s.size){
        norm += abs(t[i] - s[i])
    }
    return norm
}

public func linfDist(s: Array<Float64>, t: Array<Float64>): Float64 {
    var norm: Float64 = 0.0
    if(s.isEmpty()){
        return 0.0
    }
    norm = abs(t[0] - s[0])
    for (i in 1..s.size){
        var absDiff: Float64 = abs(t[i] - s[i]) 
        if(absDiff > norm || norm.isNaN()){
            norm = absDiff
        }
    }
    return norm
}

public func sum(x: Array<Float64>): Float64 {
    var sum: Float64 = 0.0
    for (i in 0..x.size){
        sum += x[i]
    }
    return sum
}
