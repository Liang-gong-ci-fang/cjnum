package cjnum.num.quat

import std.math.*

// Abs returns the absolute value (also called the modulus) of q.
public func nqAbs(q: NqNumber): Float64 {
    // Special cases.
    if (nqIsInf(q)) {
        return Float64.Inf
    } else if (nqIsNaN(q)) {
        return Float64.NaN
    }

    var r = q.real
    var i = q.imag
    var j = q.jmag
    var k = q.kmag
    if (r < 0.0) {
        r = -r
    }
    if (i < 0.0) {
        i = -i
    }
    if (j < 0.0) {
        j = -j
    }
    if (k < 0.0) {
        k = -k
    }
    if (r < i) {
        let temp = i
        i = r
        r = temp
    }
    if (r < j) {
        let temp = j
        j = r
        r = temp
    }
    if (r < k) {
        let temp = k
        k = r
        r = temp
    }
    if (r == 0.0) {
        return 0.0
    }
    i /= r
    j /= r
    k /= r
    return r * sqrt(1.0 + i * i + j * j + k * k)
}
