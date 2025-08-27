package cjnum.num.hyperdual

// NhNumber is a float64 precision hyperdual number.
public class NhNumber {
    let real: Float64
    let e1mag: Float64
    let e2mag: Float64
    let e1E2mag: Float64
    public init(real: Float64, e1mag: Float64, e2mag: Float64, e1E2mag: Float64) {
        this.real = real
        this.e1mag = e1mag
        this.e2mag = e2mag
        this.e1E2mag = e1E2mag
    }
}

public let negZero = Float64.fromBits(1 << 63)

// TODO func (d Number) Format(fs fmt.State, c rune)

// TODO func fmtString(fs fmt.State, c rune, prec, width int, wantPlus bool) string {

// Add returns the sum of x and y.
public func nhAdd(x: NhNumber, y: NhNumber): NhNumber {
    return NhNumber(
        x.real + y.real,
        x.e1mag + y.e1mag,
        x.e2mag + y.e2mag,
        x.e1E2mag + y.e1E2mag,
    )
}

// Sub returns the difference of x and y, x-y.
public func nhSub(x: NhNumber, y: NhNumber): NhNumber {
    return NhNumber(
        x.real - y.real,
        x.e1mag - y.e1mag,
        x.e2mag - y.e2mag,
        x.e1E2mag - y.e1E2mag,
    )
}

// Mul returns the hyperdual product of x and y.
public func nhMul(x: NhNumber, y: NhNumber): NhNumber {
    return NhNumber(
        x.real * y.real,
        x.real * y.e1mag + x.e1mag * y.real,
        x.real * y.e2mag + x.e2mag * y.real,
        x.real * y.e1E2mag + x.e1mag * y.e2mag + x.e2mag * y.e1mag + x.e1E2mag * y.real,
    )
}

// Inv returns the hyperdual inverse of d.
//
// Special cases are:
//
//	Inv(±Inf) = ±0-0ϵ₁-0ϵ₂±0ϵ₁ϵ₂
//	Inv(±0) = ±Inf-Infϵ₁-Infϵ₂±Infϵ₁ϵ₂
public func nhInv(d: NhNumber): NhNumber {
    if (d.real == 0.0) {
        return NhNumber(
            1.0 / d.real,
            -Float64.Inf,
            -Float64.Inf,
            1.0 / d.real, // Return a signed inf from a signed zero.
        )
    }
    let d2 = d.real * d.real
    return NhNumber(
        1.0 / d.real,
        -d.e1mag / d2,
        -d.e2mag / d2,
        -d.e1E2mag / d2 + 2.0 * d.e1mag * d.e2mag / (d2 * d.real),
    )
}

/*
以下TODO部分是因为原项目没有对相应方法的测试
*/
// TODO func Scale(f float64, d Number) Number {

// TODO func Abs(d Number) Number {