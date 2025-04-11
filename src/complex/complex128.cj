package cjnum.complex

import std.math.*

public class Complex128 {
    public static let NaN: Complex128 = Complex128(Float64.NaN, Float64.NaN)

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

    public operator func !=(right: Complex128) {
        return this.real != right.real || this.imag != right.imag
    }
    
    public operator func -(right: Complex128) {
        return Complex128((this.real - right.real), (this.imag - right.imag))
    }
    // TODO: add more operators and functions
}
