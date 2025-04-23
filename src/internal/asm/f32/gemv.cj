package cjnum.internal.asm.f32

public func germN(m: Int64, n: Int64, alpha: Float32, a: Array<Float32>, lda: Int64, x: Array<Float32>, incX: Int64, beta: Float32, y: Array<Float32>, incY: Int64) {
    var kx: Int64 = 0
    var ky: Int64 = 0
    if(incX < 0) {
        kx = -(n - 1) * incX
    }
    if (incY < 0) {
        ky = -(m - 1) * incY
    }

    if (incX == 1 && incY == 1) {
        if (beta == 0.0) {
            for (i in 0..m) {
                y[i] = alpha * dotUnitary(a[lda * i..lda * i + n], x)
            }
            return
        }
        for (i in 0..m) {
            y[i] = y[i] * beta + alpha * dotUnitary(a[lda * i..lda * i + n], x)
        }
        return
    }
    var iy = ky
    if(beta == 0.0){
        for (i in 0..m) {
            y[iy] = alpha * dotInc(x, a[lda * i..lda * i + n], n, incX, 1, kx, 0)
            iy += incY
        }
        return
    }
    for (i in 0..m) {
        y[iy] = y[iy] * beta + alpha * dotInc(x, a[lda * i..lda * i + n], n, incX, 1, kx, 0)
        iy += incY
    }
}

public func gemvT(m: Int64, n: Int64, alpha:Float32, a: Array<Float32>, lda: Int64, x: Array<Float32>, incX: Int64, beta: Float32, y: Array<Float32>, incY: Int64) {
    var kx: Int64 = 0
    var ky: Int64 = 0
    if (incX < 0) {
        kx = -(m - 1) * incX
    }
    if (incY < 0) {
        ky = -(n - 1) * incY
    }

    if (beta == 0.0) {
        if (incY == 1) {
            for (i in 0..y.size) {
                y[i] =0.0
            }
        } else {
            var iy = ky
            for (i in 0..n) {
                y[iy] = 0.0
                iy += incY
            }
        }
    } else if (incY < 0) {
        scalInc(beta, y, n, -incY)
    } else if (incY == 1) {
        scalUnitary(beta, y[..n])
    } else {
        scalInc(beta, y, n, incY)
    }

    if (incX == 1 && incY == 1) {
        for (i in 0..m) {
            axpyUnitaryTo(y, alpha*x[i], a[lda * i..lda * i + n], y)
        }
        return
    }

    var ix = kx
    for (i in 0..m) {
        axpyInc(alpha * x[ix], a[lda * i..lda * i + n], y, n, 1, incY, 0, ky)
        ix += incX
    }
}