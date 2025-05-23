package cjnum.num.dual

import std.math.*

// TODO public interface State {

public class Number {
    public let real: Float64
    public let emag: Float64

    public init(real: Float64, emag: Float64) {
        this.real = real
        this.emag = emag
    }

    // TODO func (d Number) Format(fs fmt.State, c rune) {
}

// TODO func fmtString(fs fmt.State, c rune, prec, width int, wantPlus bool) string

// Add returns the sum of x and y.
public func ndAdd(x: Number, y: Number): Number {
    return Number(
        x.real + y.real,
        x.emag + y.emag,
    )
}

// Mul returns the dual product of x and y.
public func ndMul(x: Number, y: Number): Number {
    return Number(
        x.real * y.real,
        x.real * y.emag + x.emag * y.real,
    )
}

// Inv returns the dual inverse of d.
//
// Special cases are:
//
//	Inv(±Inf) = ±0-0ϵ
//	Inv(±0) = ±Inf-Infϵ
public func ndInv(d: Number): Number {
    let d2 = d.real * d.real
    return Number(
        1.0 / d.real,
        -d.emag / d2,
    )
}
