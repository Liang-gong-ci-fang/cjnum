package cjnum.internal.asm.f64

public func ger(m: Int64, n: Int64, alpha: Float64, x: Array<Float64>, incX: Int64, y: Array<Float64>, incY: Int64,
    a: Array<Float64>, lda: Int64) {
    if (incX == 1 && incY == 1) {
        var x_ = x.slice(0, m)
        var y_ = y.slice(0, n)
        for (i in 0..x_.size) {
            axpyUnitary(alpha * x_[i], y_, a[lda * i..lda * i + n])
        }
        return
    }

    var ky: Int64 = 0
    var kx: Int64 = 0
    if (incY < 0) {
        ky = -(n - 1) * incY
    }
    if (incX < 0) {
        kx = -(m - 1) * incX
    }

    var ix = kx

    for (i in 0..m) {
        axpyInc(alpha * x[ix], y, a[lda * i..lda * i + n], n, incY, 1, ky, 0)
        ix += incX
    }
}

public func gemvN(m: Int64, n: Int64, alpha: Float64, a: Array<Float64>, lda: Int64, x: Array<Float64>, incX: Int64,
    beta: Float64, y: Array<Float64>, incY: Int64) {
    let kx = match {
        case incX < 0 => -(n - 1) * incX
        case _ => 0
    }
    var ky = match {
        case incY < 0 => -(m - 1) * incY
        case _ => 0
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

    if (beta == 0.0) {
        for (i in 0..m) {
            y[ky] = alpha * dotInc(x, a[lda * i..lda * i + n], n, incX, 1, kx, 0)
            ky += incY
        }
        return
    }

    for (i in 0..m) {
        y[ky] = y[ky] * beta + alpha * dotInc(x, a[lda * i..lda * i + n], n, incX, 1, kx, 0)
        ky += incY
    }
}

public func gemvT(m: Int64, n: Int64, alpha: Float64, a: Array<Float64>, lda: Int64, x: Array<Float64>, incX: Int64,
    beta: Float64, y: Array<Float64>, incY: Int64) {
    var ky: Int64 = 0
    var kx: Int64 = 0
    if (incY < 0) {
        ky = -(n - 1) * incY
    }
    if (incX < 0) {
        kx = -(m - 1) * incX
    }

    if (beta == 0.0) { // beta == 0 is special-cased to memclear
        if (incY == 1) {
            for (i in 0..y.size) {
                y[i] = 0.0
            }
        } else {
            var iy = ky
            for (_ in 0..n) {
                y[iy] = 0.0
                iy += incY
            }
        }
    } else if (incY < 0) {
        scalInc(beta, y, n, -incY)
    } else if (incY == 1) {
        scalUnitary(beta, y[0..n])
    } else {
        scalInc(beta, y, n, incY)
    }

	if (incX == 1 && incY == 1) {
		for (i in 0..m) {
            let y_ = y
			axpyUnitaryTo(y, alpha * x[i], a[lda * i..lda * i + n], y_)
		}
		return
	}

	var ix = kx
	for (i in 0..m) {
		axpyInc(alpha * x[ix], a[lda * i..lda * i + n], y, n, 1, incY, 0, ky)
		ix += incX
	}
}
