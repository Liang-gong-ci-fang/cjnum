package cjnum.num.dual

import std.math.*

// TODO public interface State {

public class NdNumber {
    public let real: Float64
    public let emag: Float64

    public init(real: Float64, emag: Float64) {
        this.real = real
        this.emag = emag
    }

    // TODO func (d NdNumber) Format(fs fmt.State, c rune) {
}

// TODO func fmtString(fs fmt.State, c rune, prec, width int, wantPlus bool) string

// Add returns the sum of x and y.
public func ndAdd(x: NdNumber, y: NdNumber): NdNumber {
    return NdNumber(
        x.real + y.real,
        x.emag + y.emag,
    )
}

// Mul returns the dual product of x and y.
public func ndMul(x: NdNumber, y: NdNumber): NdNumber {
    return NdNumber(
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
public func ndInv(d: NdNumber): NdNumber {
    let d2 = d.real * d.real
    return NdNumber(
        1.0 / d.real,
        -d.emag / d2,
    )
}

/*
以下的TODO部分因为缺乏相关的测试代码支持所以没有实现
*/
// func ndScale(f: Float64, d: NdNumber): NdNumber {
//     return NdNumber(f * d.real, f * d.emag)
// }

// func ndAbs(d: NdNumber): NdNumber {
//     if (!signbit(d.real)) {
//         return d
//     }
//     return ndScale(-1.0, d)
// }