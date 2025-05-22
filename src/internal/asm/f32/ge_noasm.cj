package cjnum.internal.asm.f32

public func ger(m: Int64, n: Int64, alpha: Float32, x:Array<Float32>, incX: Int64, y: Array<Float32>, incY: Int64, a: Array<Float32>, lda: Int64) {
    
    if (incX == 1 && incY == 1) {
        let x_ = x[..m]
        let y_ = y[..n]

        for (i in 0..x_.size) {
            axpyUnitary(alpha * x_[i], y_, a[i * lda..i * lda + n])
        }
        return
    }

    var kx: Int64 = 0
    var ky: Int64 = 0
    if (incY < 0) {
        ky = -(n - 1) * incY
    }

    if (incX < 0) {
        kx = -(m - 1) * incX
    }

    var ix = kx
    for (i in 0..m) {
        axpyInc(alpha * x[ix], y, a[i * lda..i * lda + n], n, incY, 1, ky, 0)
        ix += incX
    }
}