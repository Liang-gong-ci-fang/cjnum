package cjnum.internal.asm.c128

import cjnum.complex.*

public func axpyUnitary(alpha: Complex128, x: Array<Complex128>, y: Array<Complex128>) {
    for (i in 0..x.size) {
        y[i] += alpha * x[i]
    }
}

public func axpyUnitaryTo(dst: Array<Complex128>, alpha: Complex128, x: Array<Complex128>, y: Array<Complex128>) {
    for (i in 0..x.size) {
        dst[i] = alpha * x[i] + y[i]
    }
}

public func axpyInc(alpha: Complex128, x: Array<Complex128>, y: Array<Complex128>, n: Int64, incX: Int64, incY: Int64,
    ix: Int64, iy: Int64) {
    var ix_ = ix
    var iy_ = iy
    for (_ in 0..n) {
        y[iy_] += alpha * x[ix_]
        ix_ += incX
        iy_ += incY
    }
}

public func axpyIncTo(dst: Array<Complex128>, incDst: Int64, idst: Int64, alpha: Complex128, x: Array<Complex128>,
    y: Array<Complex128>, n: Int64, incX: Int64, incY: Int64, ix: Int64, iy: Int64) {
    var ix_ = ix
    var iy_ = iy
    var idst_ = idst
    for (_ in 0..n) {
        dst[idst_] = alpha * x[ix_] + y[iy_]
        ix_ += incX
        iy_ += incY
        idst_ += incDst
    }
}

public func dscalUnitary(alpha: Float64, x: Array<Complex128>) {
    for (i in 0..x.size) {
        x[i] = Complex128(x[i].real * alpha, x[i].imag * alpha)
    }
}

public func dscalInc(alpha: Complex128, x: Array<Complex128>, n: Int64, inc: Int64) {
    var ix: Int64 = 0
    for (_ in 0..n) {
        x[ix] *= alpha
        ix += inc
    }
}

public func scalInc(alpha: Complex128, x: Array<Complex128>, n: Int64, inc: Int64) {
    var ix: Int64 = 0
    for (_ in 0..n) {
        x[ix] *= alpha
        ix += inc
    }
}

public func scalUnitary(alpha: Complex128, x: Array<Complex128>) {
    for (i in 0..x.size) {
        x[i] = alpha * x[i]
    }
}

public func dotcUnitary(x: Array<Complex128>, y: Array<Complex128>): Complex128 {
    var sum: Complex128 = Complex128(0.0)
    for (i in 0..x.size) {
        sum += y[i] * x[i].conj()
    }
    return sum
}

public func dotcInc(x: Array<Complex128>, y: Array<Complex128>, n: Int64, incX: Int64, incY: Int64, ix: Int64, iy: Int64): Complex128 {
    var ix_ = ix
    var iy_ = iy
    var sum = Complex128(0.0)
    for (_ in 0..n) {
        sum += y[iy_] * x[ix_].conj()
        ix_ += incX
        iy_ += incY
    }
    return sum
}

public func dotuUnitary(x: Array<Complex128>, y: Array<Complex128>): Complex128 {
    var sum = Complex128(0.0)
    for (i in 0..x.size) {
        sum += y[i] * x[i]
    }
    return sum
}

public func dotuInc(x: Array<Complex128>, y: Array<Complex128>, n: Int64, incX: Int64, incY: Int64, ix: Int64, iy: Int64): Complex128 {
    var ix_ = ix
    var iy_ = iy
    var sum = Complex128(0.0)
    for (_ in 0..n) {
        sum += y[iy_] * x[ix_]
        ix_ += incX
        iy_ += incY
    }
    return sum
}
