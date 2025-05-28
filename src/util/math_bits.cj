package cjnum.util

let shift: UInt64 = 64 - 11 - 1
let mask: UInt64 = 0x7FF
let bias: UInt64 = 1023

// LeadingZeros64 returns the number of leading zero bits in x; the result is 64 for x == 0.
public func leadingZeros64(x: UInt64): Int64 {
    return 64 - len64(x)
}

// Len64 returns the minimum number of bits required to represent x; the result is 0 for x == 0.
public func len64(x_: UInt64): Int64 {
    var n = 0
    var x = x_
    if (x >= 1 << 32) {
        x >>= 32
        n = 32
    }
    if (x >= 1 << 16) {
        x >>= 16
        n += 16
    }
    if (x >= 1 << 8) {
        x >>= 8
        n += 8
    }
    return n + Int64(len8tab[Int64(x)])
}

// Add64 returns the sum with carry of x, y and carry: sum = x + y + carry.
// The carry input must be 0 or 1; otherwise the behavior is undefined.
// The carryOut output is guaranteed to be 0 or 1.
//
// This function's execution time does not depend on the inputs.
public func add64(x: UInt64, y: UInt64, carry: UInt64): (UInt64, UInt64) {
    let sum = x + y + carry
    // The sum will overflow if both top bits are set (x & y) or if one of them
    // is (x | y), and a carry from the lower place happened. If such a carry
    // happens, the top bit will be 1 + 0 + 1 = 0 (&^ sum).
    let carryOut = ((x & y) | ((x | y) & !sum)) >> 63
    return (sum, carryOut)
}

// Mul64 returns the 128-bit product of x and y: (hi, lo) = x * y
// with the product bits' upper half returned in hi and the lower
// half returned in lo.
//
// This function's execution time does not depend on the inputs.
public func mul64(x: UInt64, y: UInt64): (UInt64, UInt64) {
    let mask32: UInt64 = 1 << 32 - 1
    let x0 = x & mask32
    let x1 = x >> 32
    let y0 = y & mask32
    let y1 = y >> 32
    let w0 = x0 * y0
    let t = x1 * y0 + w0 >> 32
    var w1 = t & mask32
    let w2 = t >> 32
    w1 += x0 * y1
    let hi = x1 * y1 + w2 + w1 >> 32
    let lo = x * y
    return (hi, lo)
}
