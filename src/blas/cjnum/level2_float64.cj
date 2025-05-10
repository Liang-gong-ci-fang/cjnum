package cjnum.blas.cjnum

import cjnum.internal.asm.f64.*
import cjnum.blas.*
import std.math.*

extend NFloat64Implementation {
    public func dger(m: Int64, n: Int64, alpha: Float64, x: Array<Float64>, incX: Int64, y: Array<Float64>, incY: Int64,
        a: Array<Float64>, lda: Int64): Unit {
        if (m < 0) {
            throw Exception(mLT0)
        }
        if (n < 0) {
            throw Exception(nLT0)
        }
        if (lda < max(1, n)) {
            throw Exception(badLdA)
        }
        if (incX == 0) {
            throw Exception(zeroIncX)
        }
        if (incY == 0) {
            throw Exception(zeroIncY)
        }

        if (m == 0 || n == 0) {
            return 
        }

        if ((incX > 0 && x.size <= (m - 1) * incX) || (incX < 0 && x.size <= (1 - m) * incX)) {
            throw Exception(shortX)
        }
        if ((incY > 0 && y.size <= (n - 1) * incY) || (incY < 0 && y.size <= (1 - n) * incY)) {
            throw Exception(shortY)
        }
        if (a.size < lda * (m - 1) + n) {
            throw Exception(shortA)
        }

        if (alpha == 0.0) {
            return 
        }
        
        ger(m, n, alpha, x, incX, y, incY, a, lda)
    }

    public func dgbmv(tA: Transpose, m: Int64, n: Int64, kL: Int64, kU: Int64, alpha: Float64, a: Array<Float64>, lda: Int64,
        x: Array<Float64>, incX: Int64, beta: Float64, y: Array<Float64>, incY: Int64): Unit {
        if (tA != NoTrans && tA != Trans && tA != ConjTrans) {
            throw Exception(badTranspose)
        }
        if (m < 0) {
            throw Exception(mLT0)
        }
        if (n < 0) {
            throw Exception(nLT0)
        }
        if (kL < 0) {
            throw Exception(kLLT0)
        }
        if (kU < 0) {
            throw Exception(kULT0)
        }
        if (lda < kL + kU + 1) {
            throw Exception(badLdA)
        }
        if (incX == 0) {
            throw Exception(zeroIncX)
        }
        if (incY == 0) {
            throw Exception(zeroIncY)
        }

        if (m == 0 || n == 0) {
            return
        }

        if (a.size < lda * (min(m, n + kL) - 1) + kL + kU + 1) {
            throw Exception(shortA)
        }
        var lenX = m
        var lenY = n
        if (tA == NoTrans) {
            lenX = n
            lenY = m
        }
        if ((incX > 0 && x.size <= (lenX - 1) * incX) || (incX < 0 && x.size <= (1 - lenX) * incX)) {
            throw Exception(shortX)
        }
        if ((incY > 0 && y.size <= (lenY - 1) * incY) || (incY < 0 && y.size <= (1 - lenY) * incY)) {
            throw Exception(shortY)
        }

        if (alpha == 0.0 && beta == 1.0) {
            return 
        }

        var kx = 0
        var ky = 0
        if (incX < 0) {
            kx = -(lenX - 1) * incX
        }
        if (incY < 0) {
            ky = -(lenY - 1) * incY
        }

        if (beta != 1.0) {
            if (incY == 1) {
                if (beta == 0.0) {
                    for (i in 0..lenY) {
                        y[i] = 0.0
                    }
                } else {
                    scalUnitary(beta, y[..lenY])
                }
            } else {
                var iy = ky
                if (beta == 0.0) {
                    for (i in 0..incY) {
                        y[iy] = 0.0
                        iy += incY
                    }
                } else {
                    if (incY > 0) {
                        scalInc(beta, y, lenY, incY)
                    } else {
                        scalInc(beta, y, lenY, -incY)
                    }
                }
            }
        }

        if (alpha == 0.0) {
            return 
        }

        let nCol = kU + 1 + kL

        if (tA == NoTrans) {
            var iy = ky
            if (incX == 1) {
                for (i in 0..min(m, n + kL)) {
                    let l = max(0, kL - i)
                    let u = min(nCol, n + kL - i)
                    let off = max(0, i - kL)
                    let atmp = a[i * lda + l..i * lda + u]
                    let xtmp = x[off..off + u - l]
                    var sum = 0.0
                    for (j in 0..atmp.size) {
                        sum += xtmp[j] * atmp[j]
                    }
                    y[iy] += sum * alpha
                    iy += incY
                }
                return 
            }
            for (i in 0..min(m, n + kL)) {
                let l = max(0, kL - i)
                let u = min(nCol, n + kL - i)
                let off = max(0, i - kL)
                let atmp = a[i * lda + l..i * lda + u]
                var jx = kx
                var sum = 0.0
                for (v in atmp) {
                    sum += x[off * incX + jx] * v
                    jx += incX
                }
                y[iy] += sum * alpha
                iy += incY
            }
            return 
        }

        if (incX == 1) {
            for (i in 0..min(m, n + kL)) {
                let l = max(0, kL - i)
                let u = min(nCol, n + kL - i)
                let off = max(0, i - kL)
                let atmp = a[i * lda + l..i * lda + u]
                let tmp = alpha * x[i]
                var jy = ky
                for (v in atmp) {
                    y[jy + off * incY] += tmp * v
                    jy += incY
                }
            }
            return 
        }

        var ix = kx
        for (i in 0..min(m, n + kL)) {
            let l = max(0, kL - i)
            let u = min(nCol, n + kL - i)
            let off = max(0, i - kL)
            let atmp = a[i * lda + l..i * lda + u]
            let tmp = alpha * x[ix]
            var jy = ky
            for (v in atmp) {
                y[jy + off * incY] += tmp * v
                jy += incY
            }
            ix += incX
        }
    }

    public func dgemv(tA: Transpose, m: Int64, n: Int64, alpha: Float64, a: Array<Float64>, lda: Int64, x: Array<Float64>,
        incX: Int64, beta: Float64, y: Array<Float64>, incY: Int64): Unit {
        if (tA != NoTrans && tA != Trans && tA != ConjTrans) {
            throw Exception(badTranspose)
        }
        if (m < 0) {
            throw Exception(mLT0)
        }
        if (n < 0) {
            throw Exception(nLT0)
        }
        if (lda < max(1, n)) {
            throw Exception(badLdA)
        }
        if (incX == 0) {
            throw Exception(zeroIncX)
        }
        if (incY == 0) {
            throw Exception(zeroIncY)
        }

        var lenX = m
        var lenY = n
        if (tA == NoTrans) {
            lenX = n
            lenY = m
        }

        if (m == 0 || n == 0) {
            return 
        }

        if ((incX > 0 && x.size <= (lenX - 1) * incX) || (incX < 0 && x.size <= (1 - lenX) * incX)) {
            throw Exception(shortX)
        }
        if ((incY > 0 && y.size <= (lenY - 1) * incY) || (incY < 0 && y.size <= (1 - lenY) * incY)) {
            throw Exception(shortY)
        }
        if (a.size < lda * (m - 1) + n) {
            throw Exception(shortA)
        }

        if (alpha == 0.0 && beta == 1.0) {
            return 
        }

        if (alpha == 0.0) {
            if (incY > 0) {
                nFloat64Implementation.dscal(lenY, beta, y, incY)
            } else {
                nFloat64Implementation.dscal(lenY, beta, y, -incY)
            }
            return  
        }

        if (tA == NoTrans) {
            gemvN(m, n, alpha, a, lda, x, incX, beta, y, incY)
            return 
        }
        gemvT(m, n, alpha, a, lda, x, incX, beta, y, incY)
    }
    
    public func dtrmv(ul: Uplo, tA: Transpose, d: Diag, n: Int64, a: Array<Float64>, lda: Int64, x: Array<Float64>, incX: Int64): Unit {
        if (ul != Lower && ul != Upper) {
            throw Exception(badUplo)
        }
        if (tA != NoTrans && tA != Trans && tA != ConjTrans) {
            throw Exception(badTranspose)
        }
        if (d != NonUnitDiag && d != UnitDiag) {
            throw Exception(badDiag)
        }
        if (n < 0) {
            throw Exception(nLT0)
        }
        if (lda < max(1, n)) {
            throw Exception(badLdA)
        }
        if (incX == 0) {
            throw Exception(zeroIncX)
        }

        if (n == 0) {
            return 
        }

        if (a.size < lda * (n - 1) + n) {
            throw Exception(shortA)
        }
        if ((incX > 0 && x.size <= (n - 1) * incX) || (incX < 0 && x.size <= (1 - n) * incX)) {
            throw Exception(shortX)
        }

        let nonUnit = d != UnitDiag
        if (n == 1) {
            if (nonUnit) {
                x[0] *= a[0]
            }
            return 
        }
        var kx = 0
        if (incX <= 0) {
            kx = -(n - 1) * incX
        }
        if (tA == NoTrans) {
            if (ul == Upper) {
                if (incX == 1) {
                    for (i in 0..n) {
                        let ilda = i * lda
                        var tmp = 0.0
                        if (nonUnit) {
                            tmp = a[ilda + i] * x[i]
                        } else {
                            tmp = x[i]
                        }
                        x[i] = tmp + dotUnitary(a[ilda + i + 1..ilda + n], x[i + 1..n])
                    }
                    return 
                }
                var ix = kx
                for (i in 0..n) {
                    let ilda = i * lda
                    var tmp = 0.0
                    if (nonUnit) {
                        tmp = a[ilda + i] * x[ix]
                    } else {
                        tmp = x[ix]
                    }
                    x[ix] = tmp + dotInc(x, a[ilda + i + 1..ilda + n], n - i -1, incX, 1, ix + incX, 0)
                    ix += incX
                }
                return 
            }

            if (incX == 1) {
                var i = n - 1
                while (i >= 0) {
                    let ilda = i * lda
                    var tmp = 0.0
                    if (nonUnit) {
                        tmp += a[ilda + i] * x[i]
                    } else {
                        tmp = x[i]
                    }
                    x[i] = tmp + dotUnitary(a[ilda..ilda + i], x[..i])
                    i--
                }
                return
            }
            var ix = kx + (n - 1) * incX
            var i = n - 1
            while (i >= 0) {
                let ilda = i * lda
                var tmp = 0.0
                if (nonUnit) {
                    tmp = a[ilda + i] * x[ix]
                } else {
                    tmp = x[ix]
                }
                x[ix] = tmp + dotInc(x, a[ilda.. ilda + i], i, incX, 1, kx, 0)
                ix -= incX
                i--
            }
            return 
        }

        if (ul == Upper) {
            if (incX == 1) {
                var i = n - 1
                while (i >= 0) {
                    let ilda = i * lda
                    let xi = x[i]
                    axpyUnitary(xi, a[ilda + i + 1..ilda + n], x[i + 1..n])
                    if (nonUnit) {
                        x[i] *= a[ilda + i]
                    }
                    i--
                }
                return
            }
            var ix = kx + (n - 1) * incX
            var i = n - 1
            while (i >= 0) {
                let ilda = i * lda
                let xi = x[ix]
                axpyInc(xi, a[ilda + i + 1..ilda + n], x, n - i - 1, 1, incX, 0, kx + (i + 1) * incX)
                if (nonUnit) {
                    x[ix] *= a[ilda + i]
                }
                ix -= incX
                i--
            }
            return 
        }

        if (incX == 1) {
            for (i in 0..n) {
                let ilda = i * lda
                let xi = x[i]
                axpyUnitary(xi, a[ilda..ilda + i], x[..i])
                if (nonUnit) {
                    x[i] *= a[i * lda + i]
                }
            }
            return
        }

        var ix = kx
        for (i in 0..n) {
            let ilda = i * lda
            let xi = x[ix]
            axpyInc(xi, a[ilda..ilda + i], x, i, 1, incX, 0, kx)
            if (nonUnit) {
                x[ix] *= a[ilda + i]
            }
            ix += incX
        }
    }

    public func dtrsv(ul: Uplo, tA: Transpose, d: Diag, n: Int64, a: Array<Float64>, lda: Int64, x: Array<Float64>, incX: Int64): Unit {
        if (ul != Lower && ul != Upper) {
            throw Exception(badUplo)
        }
        if (tA != NoTrans && tA != Trans && tA != ConjTrans) {
            throw Exception(badTranspose)
        }
        if (d != NonUnitDiag && d != UnitDiag) {
            throw Exception(badDiag)
        }
        if (n < 0) {
            throw Exception(nLT0)
        }
        if (lda < max(1, n)) {
            throw Exception(badLdA)
        }
        if (incX == 0) {
            throw Exception(zeroIncX)
        }

        if (n == 0) {
            return
        }

        if (a.size < lda * (n - 1) + n) {
            throw Exception(shortA)
        }
        if ((incX > 0 && x.size <= (n - 1) * incX) || (incX < 0 && x.size <= (1 - n) * incX)) {
            throw Exception(shortX)
        }

        if (n == 1) {
            if (d == NonUnitDiag) {
                x[0] /= a[0]
            }
            return 
        }

        var kx = 0
        if (incX < 0) {
            kx = -(n - 1) * incX
        }

        let nonUnitDiag = d == NonUnitDiag
        if (tA == NoTrans) {
            if (ul == Upper) {
                if (incX == 1) {
                    var i = n - 1
                    while (i >= 0) {
                        var sum = 0.0
                        let atmp = a[i * lda + i + 1..i * lda + n]
                        for (j in 0..atmp.size) {
                            let jv = i + j + 1
                            sum += x[jv] * atmp[j]
                        }
                        x[i] -= sum
                        if (nonUnitDiag) {
                            x[i] /= a[i * lda + i]
                        }
                        i--
                    }
                    return
                }
                var ix = kx + (n - 1) * incX
                var i = n - 1
                while (i >= 0) {
                    var sum = 0.0
                    var jx = ix + incX
                    let atmp = a[i * lda + i + 1..i * lda + n]
                    for (v in atmp) {
                        sum += x[jx] * v
                        jx += incX
                    }
                    x[ix] -= sum
                    if (nonUnitDiag) {
                        x[ix] /= a[i * lda + i]
                    }
                    ix -= incX
                    i--
                }
                return 
            }

            if (incX == 1) {
                for (i in 0..n) {
                    var sum = 0.0
                    let atmp = a[i * lda.. i * lda + i]
                    for (j in 0..atmp.size) {
                        sum += x[j] * atmp[j]
                    }
                    x[i] -= sum
                    if (nonUnitDiag) {
                        x[i] /= a[i * lda + i]
                    }
                }
                return 
            }
            var ix = kx
            for (i in 0..n) {
                var jx = kx
                var sum = 0.0
                let atmp = a[i * lda..i * lda + i]
                for (v in atmp) {
                    sum += x[jx] * v
                    jx += incX
                }
                x[ix] -= sum
                if (nonUnitDiag) {
                    x[ix] /= a[i * lda + i]
                }
                ix += incX
            }
            return 
        }

        if (ul == Upper) {
            if (incX == 1) {
                for (i in 0..n) {
                    if (nonUnitDiag) {
                        x[i] /= a[i * lda + i]
                    }
                    let xi = x[i]
                    let atmp = a[i * lda + i + 1..i * lda + n]
                    for (j in 0..atmp.size) {
                        let jv = j + i + 1
                        x[jv] -= atmp[j] * xi
                    }
                }
                return 
            }
            var ix = kx
            for (i in 0..n) {
                if (nonUnitDiag) {
                    x[ix] /= a[i * lda + i]
                }
                let xi = x[ix]
                var jx = kx + (i + 1) * incX
                let atmp = a[i * lda + i + 1.. i * lda + n]
                for (v in atmp) {
                    x[jx] -= v * xi
                    jx += incX
                }
                ix += incX
            }
            return 
        }

        if (incX == 1) {
            var i = n - 1
            while (i >= 0) {
                if (nonUnitDiag) {
                    x[i] /= a[i * lda + i]
                }
                let xi = x[i]
                let atmp = a[i * lda..i * lda + i]
                for (j in 0..atmp.size) {
                    x[j] -= atmp[j] * xi
                }
                i--
            }
            return
        }
        var ix = kx + (n - 1) * incX
        var i = n - 1
        while (i >= 0) {
            if (nonUnitDiag) {
                x[ix] /= a[i * lda + i]
            }
            let xi = x[ix]
            var jx = kx
            let atmp = a[i * lda.. i * lda + i]
            for (v in atmp) {
                x[jx] -= v * xi
                jx += incX
            }
            ix -= incX
            i--
        }
    }
    
    public func dsymv(ul: Uplo, n: Int64, alpha: Float64, a: Array<Float64>, lda: Int64, x: Array<Float64>, incX: Int64,
        beta: Float64, y: Array<Float64>, incY: Int64): Unit {
        if (ul != Lower && ul != Upper) {
            throw Exception(badUplo)
        }
        if (n < 0) {
            throw Exception(nLT0)
        }
        if (lda < max(1, n)) {
            throw Exception(badLdA)
        }
        if (incX == 0) {
            throw Exception(zeroIncX)
        }
        if (incY == 0) {
            throw Exception(zeroIncY)
        }

        if (n == 0) {
            return
        }

        if (a.size < lda * (n - 1) + n) {
            throw Exception(shortA)
        }
        if ((incX > 0 && x.size <= (n - 1) * incX) || (incX < 0 && x.size <= (1 - n) * incX)) {
            throw Exception(shortX)
        }
        if ((incY > 0 && y.size <= (n - 1) * incY) || (incY < 0 && y.size <= (1 - n) * incY)) {
            throw Exception(shortY)
        }

        if (alpha == 0.0 && beta == 1.0) {
            return
        }

        var kx = 0
        var ky = 0
        if (incX < 0) {
            kx = -(n - 1) * incX
        }
        if (incY < 0) {
            ky = -(n - 1) * incY
        }

        if (beta != 1.0) {
            if (incY == 1) {
                if (beta == 0.0) {
                    for (i in 0..n) {
                        y[i] = 0.0
                    }
                } else {
                    scalUnitary(beta, y[..n])
                }
            } else {
                var iy = ky
                if (beta == 0.0) {
                    for (i in 0..n) {
                        y[iy] = 0.0
                        iy += incY
                    }
                } else {
                    if (incY > 0) {
                        scalInc(beta, y, n, incY)
                    } else {
                        scalInc(beta, y, n, -incY)
                    }
                }
            }
        }

        if (alpha == 0.0) {
            return 
        }

        if (n == 1) {
            y[0] += alpha * a[0] * x[0]
            return 
        }

        if (ul == Upper) {
            if (incX == 1) {
                var iy = ky
                for (i in 0..n) {
                    let xv = x[i] * alpha
                    var sum = x[i] * a[i * lda + i]
                    var jy = ky + (i + 1) * incY
                    let atmp = a[i * lda + i + 1..i * lda + n]
                    for (j in 0..atmp.size) {
                        let jp = j + i + 1
                        sum += x[jp] * atmp[j]
                        y[jy] += xv * atmp[j]
                        jy += incY
                    }
                    y[iy] += alpha * sum
                    iy += incY
                }
                return 
            }
            var ix = kx
            var iy = ky
            for (i in 0..n) {
                let xv = x[ix] * alpha
                var sum = x[ix] * a[i * lda + i]
                var jx = kx + (i + 1) * incX
                var jy = ky + (i + 1) * incY
                let atmp = a[i * lda + i + 1..i * lda + n]
                for (v in atmp) {
                    sum += x[jx] * v
                    y[jy] += xv * v
                    jx += incX
                    jy += incY
                }
                y[iy] += alpha * sum
                ix += incX
                iy += incY
            }
            return
        }

        if (incX == 1) {
            var iy = ky
            for (i in 0..n) {
                var jy = ky
                let xv = alpha * x[i]
                let atmp = a[i * lda..i * lda + i]
                var sum = 0.0
                for (j in 0..atmp.size) {
                    sum += x[j] * atmp[j]
                    y[jy] += xv * atmp[j]
                    jy += incY
                }
                sum += x[i] * a[i * lda + i]
                sum *= alpha
                y[iy] += sum
                iy += incY
            } 
            return 
        }

        var ix = kx
        var iy = ky
        for (i in 0..n) {
            var jx = kx
            var jy = ky
            let xv = alpha * x[ix]
            var sum = 0.0
            let atmp = a[i * lda..i * lda + i]
            for (v in atmp) {
                sum += x[jx] * v
                y[jy] += xv * v
                jx += incX
                jy += incY
            }
            sum += x[ix] * a[i * lda + i]
            sum *= alpha
            y[iy] += sum
            ix += incX
            iy += incY
        }
    }

    public func dtbmv(ul: Uplo, tA: Transpose, d: Diag, n: Int64, k: Int64, a: Array<Float64>, lda: Int64, x: Array<Float64>,
        incX: Int64): Unit {
        if (ul != Lower && ul != Upper) {
            throw Exception(badUplo)
        }
        if (tA != NoTrans && tA != Trans && tA != ConjTrans) {
            throw Exception(badTranspose)
        }
        if (d != NonUnitDiag && d != UnitDiag) {
            throw Exception(badDiag)
        }
        if (n < 0) {
            throw Exception(nLT0)
        }
        if (k < 0) {
            throw Exception(kLT0)
        }
        if (lda < k + 1) {
            throw Exception(badLdA)
        }
        if (incX == 0) {
            throw Exception(zeroIncX)
        }

        if (n == 0) {
            return
        }

        if (a.size < lda * (n - 1) + k + 1) {
            throw Exception(shortA)
        }
        if ((incX > 0 && x.size <= (n - 1) * incX) || (incX < 0 && x.size <= (1 - n) * incX)) {
            throw Exception(shortX)
        }

        var kx = 0
        if (incX < 0) {
            kx = -(n - 1) * incX
        }

        let nonUnitDiag = d != UnitDiag
        if (tA == NoTrans) {
            if (ul == Upper) {
                if (incX == 1) {
                    for (i in 0..n) {
                        let u = min(1 + k, n - i)
                        var sum = 0.0
                        let atmp = a[i * lda..]
                        let xtmp = x[i..]
                        for (j in 1..u) {
                            sum += xtmp[j] * atmp[j]
                        }
                        if (nonUnitDiag) {
                            sum += xtmp[0] * atmp[0]
                        } else {
                            sum += xtmp[0]
                        }
                        x[i] = sum
                    }
                    return 
                }
                var ix = kx
                for (i in 0..n) {
                    let u = min(1 + k, n - i)
                    var sum = 0.0
                    let atmp = a[i * lda..]
                    var jx = incX
                    for (j in 1..u) {
                        sum += x[ix + jx] * atmp[j]
                        jx += incX
                    }
                    if (nonUnitDiag) {
                        sum += x[ix] * atmp[0]
                    } else {
                        sum += x[ix]
                    }
                    x[ix] = sum
                    ix += incX
                }
                return
            }

            if (incX == 1) {
                var i = n - 1
                while (i >= 0) {
                    let l = max(0, k - i)
                    let atmp = a[i * lda..]
                    var sum = 0.0
                    for (j in l..k) {
                        sum += x[i - k + j] * atmp[j]
                    }
                    if (nonUnitDiag) {
                        sum += x[i] * atmp[k]
                    } else {
                        sum += x[i]
                    }
                    x[i] = sum
                    i--
                }
                return 
            }
            var ix = kx + (n - 1) * incX
            var i = n - 1
            while (i >= 0) {
                let l = max(0, k - i)
                let atmp = a[i * lda..]
                var sum = 0.0
                var jx = l * incX
                for (j in l..k) {
                    sum += x[ix - k * incX + jx] * atmp[j]
                    jx += incX
                }
                if (nonUnitDiag) {
                    sum += x[ix] * atmp[k]
                } else {
                    sum += x[ix]
                }
                x[ix] = sum
                ix -= incX
                i--
            }
            return 
        }

        if (ul == Upper) {
            if (incX == 1) {
                var i = n -1
                while (i >= 0) {
                    var u = k + 1
                    if (i < u) {
                        u = i + 1
                    }
                    var sum = 0.0
                    for (j in 1..u) {
                        sum += x[i - j] * a[(i - j) * lda + j]
                    }
                    if (nonUnitDiag) {
                        sum += x[i] * a[i * lda]
                    } else {
                        sum += x[i]
                    }
                    x[i] = sum
                    i--
                }
                return 
            }
            var ix = kx + (n - 1) * incX
            var i = n - 1
            while (i >= 0) {
                var u = k + 1
                if (i < u) {
                    u = i + 1 
                }
                var sum = 0.0
                var jx = incX
                for (j in 1..u) {
                    sum += x[ix - jx] * a[(i - j) * lda + j]
                    jx += incX
                }
                if (nonUnitDiag) {
                    sum += x[ix] * a[i * lda]
                } else {
                    sum += x[ix]
                }
                x[ix] = sum
                ix -= incX
                i--
            }
            return 
        }

        if (incX == 1) {
            for (i in 0..n) {
                var u = k
                if (i + k >= n) {
                    u = n - i - 1
                }
                var sum = 0.0
                for (j in 0..u) {
                    sum += x[i + j + 1] * a[(i + j + 1) * lda + k - j - 1]
                }
                if (nonUnitDiag) {
                    sum += x[i] * a[i * lda + k]
                } else {
                    sum += x[i]
                }
                x[i] = sum
            }
            return 
        }

        var ix = kx 
        for (i in 0..n) {
            var u = k
            if (i + k >= n) {
                u = n - i - 1
            }
            var sum = 0.0
            var jx = 0
            for (j in 0..u) {
                sum += x[ix + jx + incX] * a[(i + j + 1) * lda + k - j - 1]
                jx += incX
            } 
            if (nonUnitDiag) {
                sum += x[ix] * a[i * lda + k]
            } else {
                sum += x[ix]
            }
            x[ix] = sum
            ix += incX
        }
    }
    
    public func dtpmv(ul: Uplo, tA: Transpose, d: Diag, n: Int64, ap: Array<Float64>, x: Array<Float64>, incX: Int64): Unit {
        if (ul != Lower && ul != Upper) {
            throw Exception(badUplo)
        }
        if (tA != NoTrans && tA != Trans && tA != ConjTrans) {
            throw Exception(badTranspose)
        }
        if (d != NonUnitDiag && d != UnitDiag) {
            throw Exception(badDiag)
        }
        if (n < 0) {
            throw Exception(nLT0)
        }
        if (incX == 0) {
            throw Exception(zeroIncX)
        }

        if (n == 0) {
            return
        }

        if (ap.size < n * (n + 1) / 2) {
            throw Exception(shortAP)
        }

        if ((incX > 0 && x.size <= (n - 1) * incX) || (incX < 0 && x.size <= (1 - n) * incX)) {
            throw Exception(shortX)
        }

        var kx = 0
        if (incX < 0) {
            kx = -(n - 1) * incX
        }

        let nonUnitDiag = d == NonUnitDiag
        var offset = 0
        if (tA == NoTrans) {
            if (ul == Upper) {
                if (incX == 1) {
                    for (i in 0..n) {
                        var xi = x[i]
                        if (nonUnitDiag) {
                            xi *= ap[offset]
                        }
                        let atmp = ap[offset + 1..offset + n - i]
                        let xtmp = x[i + 1..]
                        for (j in 0..atmp.size) {
                            xi += atmp[j] * xtmp[j]
                        }
                        x[i] = xi
                        offset += n - i
                    }
                    return
                }

                var ix = kx
                for (i in 0..n) {
                    var xix = x[ix]
                    if (nonUnitDiag) {
                        xix *= ap[offset]
                    }
                    let atmp = ap[offset + 1..offset + n - i]
                    var jx = kx + (i + 1) * incX
                    for (v in atmp) {
                        xix += v * x[jx]
                        jx += incX
                    }
                    x[ix] = xix
                    offset += n - i
                    ix += incX
                }
                return
            }

            if (incX == 1) {
                var offset = n * (n + 1) / 2 - 1
                var i = n - 1
                while (i >= 0) {
                    var xi = x[i]
                    if (nonUnitDiag) {
                        xi *= ap[offset]
                    }
                    let atmp = ap[offset - i..offset]
                    for (j in 0..atmp.size) {
                        xi += atmp[j] * x[j]
                    }
                    x[i] = xi
                    offset -= i + 1
                    i--
                }
                return
            }
            var ix = kx + (n - 1) * incX
            var offset = n * (n + 1) / 2 - 1
            var i = n - 1
            while (i >= 0) {
                var xix = x[ix]
                if (nonUnitDiag) {
                    xix *= ap[offset]
                }
                let atmp = ap[offset - i..offset]
                var jx = kx
                for (v in atmp) {
                    xix += v * x[jx]
                    jx += incX
                }
                x[ix] = xix
                offset -= i + 1
                ix -= incX
                i--
            }
            return 
        }

        if (ul == Upper) {
            if (incX == 1) {
                var offset = n * (n + 1) / 2 - 1
                var i = n - 1
                while (i >= 0) {
                    var xi = x[i]
                    let atmp = ap[offset + 1..offset + n - i]
                    let xtmp = x[i + 1..]
                    for (j in 0..atmp.size) {
                        xtmp[j] += atmp[j] * xi
                    }
                    if (nonUnitDiag) {
                        x[i] *= ap[offset]
                    }
                    offset -= n - i + 1
                    i--
                }
                return
            }
            var ix = kx + (n - 1) * incX
            var offset = n * (n + 1) / 2 - 1
            var i = n - 1
            while (i >= 0) {
                var xix = x[ix]
                var jx = kx + (i + 1) * incX
                let atmp = ap[offset + 1..offset + n - i]
                for (v in atmp) {
                    x[jx] += v * xix
                    jx += incX
                }
                if (nonUnitDiag) {
                    x[ix] *= ap[offset]
                }
                offset -= n - i + 1
                ix -= incX
                i--
            }
            return 
        }

        if (incX == 1) {
            for (i in 0..n) {
                var xi = x[i]
                let atmp = ap[offset - i..offset]
                for (j in 0..atmp.size) {
                    x[j] += atmp[j] * xi
                }
                if (nonUnitDiag) {
                    x[i] *= ap[offset]
                }
                offset += i + 2
            }
            return 
        }

        var ix = kx
        for (i in 0..n) {
            var xix = x[ix]
            var jx = kx
            let atmp = ap[offset - i..offset]
            for (v in atmp) {
                x[jx] += v * xix
                jx += incX
            }
            if (nonUnitDiag) {
                x[ix] *= ap[offset]
            }
            ix += incX
            offset += i + 2
        }
    }
    
    public func dtbsv(ul: Uplo, tA: Transpose, d: Diag, n: Int64, k: Int64, a: Array<Float64>, lda: Int64, x: Array<Float64>,
        incX: Int64): Unit {
        if (ul != Lower && ul != Upper) {
            throw Exception(badUplo)
        }
        if (tA != NoTrans && tA != Trans && tA != ConjTrans) {
            throw Exception(badTranspose)
        }
        if (d != NonUnitDiag && d != UnitDiag) {
            throw Exception(badDiag)
        }
        if (n < 0) {
            throw Exception(nLT0)
        }
        if (k < 0) {
            throw Exception(kLT0)
        }
        if (lda < k + 1) {
            throw Exception(badLdA)
        }
        if (incX == 0) {
            throw Exception(zeroIncX)
        }

        if (n == 0) {
            return
        }

        if (a.size < lda * (n - 1) + k + 1) {
            throw Exception(shortA)
        }

        if ((incX > 0 && x.size <= (n - 1) * incX) || (incX < 0 && x.size <= (1 - n) * incX)) {
            throw Exception(shortX)
        }

        var kx = 0
        if (incX < 0) {
            kx = -(n - 1) * incX
        }

        let nonUnitDiag = d == NonUnitDiag
        if (tA == NoTrans) {
            if (ul == Upper) {
                if (incX == 1) {
                    var i = n - 1
                    while (i >= 0) {
                        var bands = k
                        if (i + bands >= n) {
                            bands = n - i - 1
                        }
                        let atmp = a[i * lda + 1..]
                        let xtmp = x[i + 1..i + bands + 1]
                        var sum = 0.0
                        for (j in 0..xtmp.size) {
                            sum += xtmp[j] * atmp[j]
                        }
                        x[i] -= sum
                        if (nonUnitDiag) {
                            x[i] /= a[i * lda]
                        }
                        i--
                    }
                    return 
                }

                var ix = kx + (n - 1) * incX
                var i = n - 1
                while (i >= 0) {
                    var max_ = k + 1
                    if (i + max_ > n) {
                        max_ = n - i
                    }
                    let atmp = a[i * lda..]
                    var jx = 0
                    var sum = 0.0
                    for (j in 1..max_) {
                        jx += incX
                        sum += x[ix + jx] * atmp[j]
                    }
                    x[ix] -= sum
                    if (nonUnitDiag) {
                        x[ix] /= atmp[0]
                    }
                    ix -= incX
                    i--
                }
                return 
            }

            if (incX == 1) {
                for (i in 0..n) {
                    var bands = k
                    if (i - k < 0) {
                        bands = i
                    }
                    let atmp = a[i * lda + k - bands..]
                    let xtmp = x[i - bands..i]
                    var sum = 0.0
                    for (j in 0..xtmp.size) {
                        sum += xtmp[j] * atmp[j]
                    }
                    x[i] -= sum
                    if (nonUnitDiag) {
                        x[i] /= atmp[bands]
                    }
                }
                return 
            }

            var ix = kx
            for (i in 0..n) {
                var bands = k
                if (i - k < 0) {
                    bands = i
                } 
                let atmp = a[i * lda + k - bands..]
                var sum = 0.0
                var jx = 0
                for (j in 0..bands) {
                    sum += x[ix - bands * incX + jx] * atmp[j]
                    jx += incX
                }
                x[ix] -= sum
                if (nonUnitDiag) {
                    x[ix] /= atmp[bands]
                }
                ix += incX
            }
            return 
        }

        if (ul == Upper) {
            if (incX == 1) {
                for (i in 0..n) {
                    var bands = k
                    if (i - k < 0) {
                        bands = i
                    }
                    var sum = 0.0
                    for (j in 0..bands) {
                        sum += x[i - bands + j] * a[(i - bands + j) * lda + bands - j]
                    }
                    x[i] -= sum
                    if (nonUnitDiag) {
                        x[i] /= a[i * lda]
                    }
                }
                return 
            }

            var ix = kx
            for (i in 0..n) {
                var bands = k
                if (i - k < 0) {
                    bands = i
                }
                var sum = 0.0
                var jx = 0
                for (j in 0..bands) {
                    sum += x[ix - bands * incX + jx] * a[(i - bands + j) * lda + bands - j]
                    jx += incX
                }
                x[ix] -= sum
                if (nonUnitDiag) {
                    x[ix] /= a[i * lda]
                }
                ix += incX
            }
            return 
        }

        if (incX == 1) {
            var i = n - 1
            while (i >= 0) {
                var bands = k
                if (i + bands >= n) {
                    bands = n - i - 1
                }
                var sum = 0.0
                let xtmp = x[i + 1..i + 1 + bands]
                for (j in 0..xtmp.size) {
                    sum += xtmp[j] * a[(i + j + 1) * lda + k - j - 1]
                }
                x[i] -= sum
                if (nonUnitDiag) {
                    x[i] /= a[i * lda + k]
                }
                i--
            }
            return 
        }

        var ix = kx + (n - 1) * incX
        var i = n - 1
        while (i >= 0) {
            var bands = k
            if (i + bands >= n) {
                bands = n - i - 1
            }
            var sum = 0.0
            var jx = 0
            for (j in 0..bands) {
                sum += x[ix + jx + incX] * a[(i + j + 1) * lda + k - j - 1]
                jx += incX
            }
            x[ix] -= sum
            if (nonUnitDiag) {
                x[ix] /= a[i * lda + k]
            }
            ix -= incX
            i--
        }
    }
    
    public func dsbmv(ul: Uplo, n: Int64, k: Int64, alpha: Float64, a: Array<Float64>, lda: Int64, x: Array<Float64>,
        incX: Int64, beta: Float64, y: Array<Float64>, incY: Int64): Unit {
        if (ul != Lower && ul != Upper) {
            throw Exception(badUplo)
        }
        if (n < 0) {
            throw Exception(nLT0)
        }
        if (k < 0) {
            throw Exception(kLT0)
        }
        if (lda < k + 1) {
            throw Exception(badLdA)
        }
        if (incX == 0) {
            throw Exception(zeroIncX)
        }
        if (incY == 0) {
            throw Exception(zeroIncY)
        }

        if (n == 0) {
            return
        }

        if (a.size < lda * (n - 1) + k + 1) {
            throw Exception(shortA)
        }
        if ((incX > 0 && x.size <= (n - 1) * incX) || (incX < 0 && x.size <= (1 - n) * incX)) {
            throw Exception(shortX)
        }
        if ((incY > 0 && y.size <= (n - 1) * incY) || (incY < 0 && y.size <= (1 - n) * incY)) {
            throw Exception(shortY)
        }

        if (alpha == 0.0 && beta == 1.0) {
            return
        }

        let lenX = n
        let lenY = n
        var kx = 0
        var ky = 0
        if (incX < 0) {
            kx = -(lenX - 1) * incX
        }
        if (incY < 0) {
            ky = -(lenY - 1) * incY
        }

        if (beta != 1.0) {
            if (incY == 1) {
                if (beta == 0.0) {
                    for (i in 0..n) {
                        y[i] = 0.0
                    }
                } else {
                    scalUnitary(beta, y[..n])
                }
            } else {
                var iy = ky
                if (beta == 0.0) { 
                    for (i in 0..n) {
                        y[iy] = 0.0
                        iy += incY
                    }
                } else {
                    if (incY > 0) {
                        scalInc(beta, y, n, incY)
                    } else {
                        scalInc(beta, y, n, -incY)
                    }
                }
            }
        }

        if (alpha == 0.0) {
            return 
        }

        if (ul == Upper) {
            if (incX == 1) {
                var iy = ky
                for (i in 0..n) {
                    let atmp = a[i * lda..]
                    let tmp = alpha * x[i]
                    var sum = tmp * atmp[0]
                    let u = min(k, n - i - 1)
                    var jy = incY
                    for (j in 1..=u) {
                        let v = atmp[j]
                        sum += alpha * x[i + j] * v
                        y[iy + jy] += tmp * v
                        jy += incY
                    }
                    y[iy] += sum
                    iy += incY
                }
                return 
            }
            var ix = kx
            var iy = ky
            for (i in 0..n) {
                let atmp = a[i * lda..]
                let tmp = alpha * x[ix]
                var sum = tmp * atmp[0]
                let u = min(k, n - i - 1)
                var jx = incX
                var jy = incY
                for (j in 1..=u) {
                    let v = atmp[j]
                    sum += alpha * x[ix + jx] * v
                    y[iy + jy] += tmp * v
                    jx += incX
                    jy += incY
                }
                y[iy] += sum
                ix += incX
                iy += incY
            }
            return
        }

        if (incX == 1) {
            var iy = ky
            for (i in 0..n) {
                let l = max(0, k - i)
                let tmp = alpha * x[i]
                var jy = l * incY
                let atmp = a[i * lda..]
                for (j in l..k) {
                    let v = atmp[j]
                    y[iy] += alpha * v * x[i - k + j]
                    y[iy - k * incY + jy] += tmp * v
                    jy += incY
                }
                y[iy] += tmp * atmp[k]
                iy += incY
            }
            return
        }

        var ix = kx
        var iy = ky
        for (i in 0..n) {
            let l = max(0, k - i)
            let tmp = alpha * x[ix]
            var jx = l * incX
            var jy = l * incY
            let atmp = a[i * lda..]
            for (j in l..k) {
                let v = atmp[j]
                y[iy] += alpha * v * x[ix - k * incX + jx]
                y[iy - k * incY + jy] += tmp * v
                jx += incX
                jy += incY
            }
            y[iy] += tmp * atmp[k]
            ix += incX
            iy += incY
        }
    }

    public func dsyr(ul: Uplo, n: Int64, alpha: Float64, x: Array<Float64>, incX: Int64, a: Array<Float64>, lda: Int64): Unit {
        if (ul != Lower && ul != Upper) {
            throw Exception(badUplo)
        }
        if (n < 0) {
            throw Exception(nLT0)
        }
        if (lda < max(1, n)) {
            throw Exception(badLdA)
        }
        if (incX == 0) {
            throw Exception(zeroIncX)
        }

        if (n == 0) {
            return
        }

        if ((incX > 0 && x.size <= (n - 1) * incX) || (incX < 0 && x.size <= (1 - n) * incX)) {
            throw Exception(shortX)
        }

        if (a.size < lda * (n - 1) + n) {
            throw Exception(shortA)
        }

        if (alpha == 0.0) {
            return
        }

        let lenX = n
        var kx = 0
        if (incX < 0) {
            kx = -(lenX - 1) * incX
        }

        if (ul == Upper) {
            if (incX == 1) {
                for (i in 0..n) {
                    let tmp = x[i] * alpha
                    if (tmp != 0.0) {
                        let atmp = a[i * lda + i..i * lda + n]
                        let xtmp = x[i..n]
                        for (j in 0..xtmp.size) {
                            atmp[j] += xtmp[j] * tmp
                        }
                    }
                }
                return 
            }

            var ix = kx
            for (i in 0..n) {
                let tmp = x[ix] * alpha
                if (tmp != 0.0) {
                    var jx = ix
                    let atmp = a[i * lda..]
                    for (j in i..n) {
                        atmp[j] += x[jx] * tmp
                        jx += incX
                    }
                }
                ix += incX
            }
            return 
        }

        if (incX == 1) {
            for (i in 0..n) {
                let tmp = x[i] * alpha
                if (tmp != 0.0) {
                    let atmp = a[i * lda..]
                    let xtmp = x[..i + 1]
                    for (j in 0..xtmp.size) {
                        atmp[j] += tmp * xtmp[j]
                    }
                }
            }
            return
        }

        var ix = kx
        for (i in 0..n) {
            let tmp = x[ix] * alpha
            if (tmp != 0.0) {
                let atmp = a[i * lda..]
                var jx = kx
                for (j in 0..i + 1) {
                    atmp[j] += tmp * x[jx]
                    jx += incX
                }
            }
            ix += incX
        }
    }

    public func dsyr2(ul: Uplo, n: Int64, alpha: Float64, x: Array<Float64>, incX: Int64, y: Array<Float64>, incY: Int64,
        a: Array<Float64>, lda: Int64): Unit {
        if (ul != Lower && ul != Upper) {
            throw Exception(badUplo)
        }
        if (n < 0) {
            throw Exception(nLT0)
        }
        if (lda < max(1, n)) {
            throw Exception(badLdA)
        }
        if (incX == 0) {
            throw Exception(zeroIncX)
        }
        if (incY == 0) {
            throw Exception(zeroIncY)
        }

        if (n == 0) {
            return
        }

        if ((incX > 0 && x.size <= (n - 1) * incX) || (incX < 0 && x.size <= (1 - n) * incX)) {
            throw Exception(shortX)
        }
        if ((incY > 0 && y.size <= (n - 1) * incY) || (incY < 0 && y.size <= (1 - n) * incY)) {
            throw Exception(shortY)
        }

        if (a.size < lda * (n - 1) + n) {
            throw Exception(shortA)
        }

        if (alpha == 0.0) {
            return
        }

        var kx = 0
        var ky = 0
        if (incX < 0) {
            kx = -(n - 1) * incX
        }
        if (incY < 0) {
            ky = -(n - 1) * incY
        }

        if (ul == Upper) {
            if (incX == 1 && incY == 1) {
                for (i in 0..n) {
                    let xi = x[i]
                    let yi = y[i]
                    let atmp = a[i * lda..]
                    for (j in i..n) {
                        atmp[j] += alpha * (xi * y[j] + x[j] * yi)
                    }
                }
                return
            }

            var ix = kx
            var iy = ky
            for (i in 0..n) {
                var jx = kx + i * incX
                var jy = ky + i * incY
                let xi = x[ix]
                let yi = y[iy]
                let atmp = a[i * lda..]
                for (j in i..n) {
                    atmp[j] += alpha * (xi * y[jy] + x[jx] * yi)
                    jx += incX
                    jy += incY
                }
                ix += incX
                iy += incY
            }
            return 
        }
        
        if (incX == 1 && incY == 1) {
            for (i in 0..n) {
                let xi = x[i]
                let yi = y[i]
                let atmp = a[i * lda..]
                for (j in 0..=i) {
                    atmp[j] += alpha * (xi * y[j] + x[j] * yi)
                }
            }
            return
        }

        var ix = kx
        var iy = ky
        for (i in 0..n) {
            var jx = kx
            var jy = ky
            let xi = x[ix]
            let yi = y[iy]
            let atmp = a[i * lda..]
            for (j in 0..=i) {
                atmp[j] += alpha * (xi * y[jy] + x[jx] * yi)
                jx += incX
                jy += incY
            }
            ix += incX
            iy += incY
        }
    }

    public func dtpsv(ul: Uplo, tA: Transpose, d: Diag, n: Int64, ap: Array<Float64>, x: Array<Float64>, incX: Int64): Unit {
        if (ul != Lower && ul != Upper) {
            throw Exception(badUplo)
        }
        if (tA != NoTrans && tA != Trans && tA != ConjTrans) {
            throw Exception(badTranspose)
        }
        if (d != NonUnitDiag && d != UnitDiag) {
            throw Exception(badDiag)
        }
        if (n < 0) {
            throw Exception(nLT0)
        }
        if (incX == 0) {
            throw Exception(zeroIncX)
        }

        if (n == 0) {
            return
        }

        if (ap.size < n * (n + 1) / 2) {
            throw Exception(shortAP)
        }

        if ((incX > 0 && x.size <= (n - 1) * incX) || (incX < 0 && x.size <= (1 - n) * incX)) {
            throw Exception(shortX)
        }

        var kx = 0
        if (incX < 0) {
            kx = -(n - 1) * incX
        }

        let nonUnitDiag = d == NonUnitDiag
        var offset = 0
        if (tA == NoTrans) {
            if (ul == Upper) {
                offset = n * (n + 1) / 2 - 1
                if (incX == 1) {
                    var i = n - 1
                    while (i >= 0) {
                        let atmp = ap[offset + 1..offset + n - i]
                        let xtmp = x[i + 1..]
                        var sum = 0.0
                        for (j in 0..atmp.size) {
                            sum += atmp[j] * xtmp[j]
                        }
                        x[i] -= sum
                        if (nonUnitDiag) {
                            x[i] /= ap[offset]
                        }
                        offset -= n - i + 1
                        i--
                    }
                    return
                }

                var ix = kx + (n - 1) * incX
                var i = n - 1
                while (i >= 0) {
                    let atmp = ap[offset + 1..offset + n - i]
                    var jx = kx + (i + 1) * incX
                    var sum = 0.0
                    for (v in atmp) {
                        sum += v * x[jx]
                        jx += incX
                    }
                    x[ix] -= sum
                    if (nonUnitDiag) {
                        x[ix] /= ap[offset]
                    }
                    ix -= incX
                    offset -= n - i + 1
                    i--
                }
                return
            }

            if (incX == 1) {
                for (i in 0..n) {
                    let atmp = ap[offset - i..offset]
                    var sum = 0.0
                    for (j in 0..atmp.size) {
                        sum += atmp[j] * x[j]
                    }
                    x[i] -= sum
                    if (nonUnitDiag) {
                        x[i] /= ap[offset]
                    }
                    offset += i + 2
                }
                return
            }

            var ix = kx
            for (i in 0..n) {
                var jx = kx
                let atmp = ap[offset - i..offset]
                var sum = 0.0
                for (v in atmp) {
                    sum += v * x[jx]
                    jx += incX
                }
                x[ix] -= sum
                if (nonUnitDiag) {
                    x[ix] /= ap[offset]
                }
                ix += incX
                offset += i + 2
            }
            return 
        }

        if (ul == Upper) {
            if (incX == 1) {
                for (i in 0..n) {
                    if (nonUnitDiag) {
                        x[i] /= ap[offset]
                    }
                    let xi = x[i]
                    let atmp = ap[offset + 1..offset + n - i]
                    let xtmp = x[i + 1..]
                    for (j in 0..atmp.size) {
                        xtmp[j] -= atmp[j] * xi
                    }
                    offset += n - i
                }
                return 
            }

            var ix = kx
            for (i in 0..n) {
                if (nonUnitDiag) {
                    x[ix] /= ap[offset]
                }
                let xix = x[ix]
                let atmp = ap[offset + 1..offset + n - i]
                var jx = kx + (i + 1) * incX
                for (v in atmp) {
                    x[jx] -= v * xix
                    jx += incX
                }
                ix += incX
                offset += n - i
            }
            return 
        }

        if (incX == 1) {
            offset = n * (n + 1) / 2 - 1
            var i = n - 1
            while (i >= 0) {
                if (nonUnitDiag) {
                    x[i] /= ap[offset]
                }
                let xi = x[i]
                let atmp = ap[offset - i..offset]
                for (j in 0..atmp.size) {
                    x[j] -= atmp[j] * xi
                }
                offset -= i + 1
                i--
            }
            return 
        }

        var ix = kx + (n - 1) * incX
        offset = n * (n + 1) / 2 - 1
        var i = n - 1
        while (i >= 0) {
            if (nonUnitDiag) {
                x[ix] /= ap[offset]
            }
            let xix = x[ix]
            let atmp = ap[offset - i..offset]
            var jx = kx
            for (v in atmp) {
                x[jx] -= v * xix
                jx += incX
            }
            ix -= incX
            offset -= i + 1
            i--
        }
    }
    
    public func dspmv(ul: Uplo, n: Int64, alpha: Float64, ap: Array<Float64>, x: Array<Float64>, incX: Int64, beta: Float64,
        y: Array<Float64>, incY: Int64): Unit {
        if (ul != Lower && ul != Upper) {
            throw Exception(badUplo)
        }
        if (n < 0) {
            throw Exception(nLT0)
        }
        if (incX == 0) {
            throw Exception(zeroIncX)
        }
        if (incY == 0) {
            throw Exception(zeroIncY)
        }

        if (n == 0) {
            return
        }

        if (ap.size < n * (n + 1) / 2) {
            throw Exception(shortAP)
        }
        if ((incX > 0 && x.size <= (n - 1) * incX) || (incX < 0 && x.size <= (1 - n) * incX)) {
            throw Exception(shortX)
        }
        if ((incY > 0 && y.size <= (n - 1) * incY) || (incY < 0 && y.size <= (1 - n) * incY)) {
            throw Exception(shortY)
        }

        if (alpha == 0.0 && beta == 1.0) {
            return
        }

        var kx = 0
        var ky = 0
        if (incX < 0) {
            kx = -(n - 1) * incX
        }
        if (incY < 0) {
            ky = -(n - 1) * incY
        }

        if (beta != 1.0) {
            if (incY == 1) {
                if (beta == 0.0) {
                    for (i in 0..n) {
                        y[i] = 0.0
                    }
                } else {
                    scalUnitary(beta, y[..n])
                }
            } else {
                var iy = ky
                if (beta == 0.0) {
                    for (i in 0..n) {
                        y[iy] = 0.0
                        iy += incY
                    }
                } else {
                    if (incY > 0) {
                        scalInc(beta, y, n, incY)
                    } else {
                        scalInc(beta, y, n, -incY)
                    }
                }
            }
        }

        if (alpha == 0.0) {
            return
        }

        if (n == 1) {
            y[0] += alpha * ap[0] * x[0]
            return
        }

        var offset = 0
        if (ul == Upper) {
            if (incX == 1) {
                var iy = ky
                for (i in 0..n) {
                    let xv = x[i] * alpha
                    var sum = ap[offset] * x[i]
                    let atmp = ap[offset + 1..offset + n - i]
                    let xtmp = x[i + 1..]
                    var jy = ky + (i + 1) * incY
                    for (j in 0..atmp.size) {
                        sum += atmp[j] * xtmp[j]
                        y[jy] += atmp[j] * xv
                        jy += incY
                    }
                    y[iy] += alpha * sum
                    iy += incY
                    offset += n - i
                }
                return 
            }

            var ix = kx
            var iy = ky
            for (i in 0..n) {
                let xv = x[ix] *alpha
                var sum = ap[offset] * x[ix]
                let atmp = ap[offset + 1..offset + n - i]
                var jx = kx + (i + 1) * incX
                var jy = ky + (i + 1) * incY
                for (v in atmp) {
                    sum += v * x[jx]
                    y[jy] += v * xv
                    jx += incX
                    jy += incY
                }
                y[iy] += alpha * sum
                ix += incX
                iy += incY
                offset += n - i
            }
            return 
        }

        if (incX == 1) {
            var iy = ky
            for (i in 0..n) {
                let xv = x[i] *alpha
                let atmp = ap[offset - i..offset]
                var jy = ky
                var sum = 0.0
                for (j in 0..atmp.size) {
                    sum += atmp[j] * x[j]
                    y[jy] += atmp[j] * xv
                    jy += incY
                }
                sum += ap[offset] * x[i]
                y[iy] += alpha * sum
                iy += incY
                offset += i + 2
            }
            return 
        }

        var ix = kx
        var iy = ky
        for (i in 0..n) {
            let xv = x[ix] * alpha
            let atmp = ap[offset - i..offset]
            var jx = kx
            var jy = ky
            var sum = 0.0
            for (v in atmp) {
                sum += v * x[jx]
                y[jy] += v * xv
                jx += incX
                jy += incY
            }
            sum += ap[offset] * x[ix]
            y[iy] += alpha * sum
            ix += incX
            iy += incY
            offset += i + 2
        }
    }
    
    public func dspr(ul: Uplo, n: Int64, alpha: Float64, x: Array<Float64>, incX: Int64, ap: Array<Float64>): Unit {
        if (ul != Lower && ul != Upper) {
            throw Exception(badUplo)
        }
        if (n < 0) {
            throw Exception(nLT0)
        }
        if (incX == 0) {
            throw Exception(zeroIncX)
        }

        if (n == 0) {
            return
        }

        if ((incX > 0 && x.size <= (n - 1) * incX) || (incX < 0 && x.size <= (1 - n) * incX)) {
            throw Exception(shortX)
        }
        if (ap.size < n * (n + 1) / 2) {
            throw Exception(shortAP)
        }

        if (alpha == 0.0) {
            return
        }

        let lenX = n
        var kx = 0
        if (incX < 0) {
            kx = -(lenX - 1) * incX
        }

        var offset = 0
        if (ul == Upper) {
            if (incX == 1) {
                for (i in 0..n) {
                    let atmp = ap[offset..]
                    let xv = alpha * x[i]
                    let xtmp = x[i..n]
                    for (j in 0..xtmp.size) {
                        atmp[j] += xv * xtmp[j]
                    }
                    offset += n - i
                }
                return 
            }

            var ix = kx
            for (i in 0..n) {
                var jx = kx + i * incX
                let atmp = ap[offset..]
                let xv = alpha * x[ix]
                for (j in 0..n - i) {
                    atmp[j] += xv * x[jx]
                    jx += incX
                }
                ix += incX
                offset += n - i
            }
            return 
        }

        if (incX == 1) {
            for (i in 0..n) {
                let atmp = ap[offset - i..]
                let xv = alpha * x[i]
                let xtmp = x[..i + 1]
                for (j in 0..xtmp.size) {
                    atmp[j] += xv * xtmp[j]
                }
                offset += i + 2
            }
            return 
        }

        var ix = kx
        for (i in 0..n) {
            var jx = kx
            let atmp = ap[offset - i..]
            let xv = alpha * x[ix]
            for (j in 0..=i) {
                atmp[j] += xv * x[jx]
                jx += incX
            }
            ix += incX
            offset += i + 2
        }
    }
    
    public func dspr2(ul: Uplo, n: Int64, alpha: Float64, x: Array<Float64>, incX: Int64, y: Array<Float64>, incY: Int64,
        ap: Array<Float64>): Unit {
        if (ul != Lower && ul != Upper) {
            throw Exception(badUplo)
        }
        if (n < 0) {
            throw Exception(nLT0)
        }
        if (incX == 0) {
            throw Exception(zeroIncX)
        }
        if (incY == 0) {
            throw Exception(zeroIncY)
        }

        if (n == 0) {
            return
        }

        if ((incX > 0 && x.size <= (n - 1) * incX) || (incX < 0 && x.size <= (1 - n) * incX)) {
            throw Exception(shortX)
        }
        if ((incY > 0 && y.size <= (n - 1) * incY) || (incY < 0 && y.size <= (1 - n) * incY)) {
            throw Exception(shortY)
        }
        if (ap.size < n * (n + 1) / 2) {
            throw Exception(shortAP)
        }

        if (alpha == 0.0) {
            return
        }

        var kx = 0
        var ky = 0
        if (incX < 0) {
            kx = -(n - 1) * incX
        }
        if (incY < 0) {
            ky = -(n - 1) * incY
        }

        var offset = 0
        if (ul == Upper) {
            if (incX == 1 && incY == 1) {
                for (i in 0..n) {
                    let atmp = ap[offset..]
                    let xi = x[i]
                    let yi = y[i]
                    let xtmp = x[i..n]
                    let ytmp = y[i..n]
                    for (j in 0..xtmp.size) {
                        atmp[j] += alpha * (xi * ytmp[j] + xtmp[j] * yi)
                    }
                    offset += n - i
                }
                return 
            }

            var ix = kx
            var iy = ky
            for (i in 0..n) {
                var jx = kx + i * incX
                var jy = ky + i * incY
                let atmp = ap[offset..]
                let xi = x[ix]
                let yi = y[iy]
                for (j in 0..n - i) {
                    atmp[j] += alpha * (xi * y[jy] + x[jx] * yi)
                    jx += incX
                    jy += incY
                }
                ix += incX
                iy += incY
                offset += n - i
            }
            return 
        }

        if (incX == 1 && incY == 1) {
            for (i in 0..n) {
                let atmp = ap[offset - i..]
                let xi = x[i]
                let yi = y[i]
                let xtmp = x[..i + 1]
                for (j in 0..xtmp.size) {
                    atmp[j] += alpha * (xi * y[j] + xtmp[j] * yi)
                }
                offset += i + 2
            }
            return
        }

        var ix = kx
        var iy = ky
        for (i in 0..n) {
            var jx = kx
            var jy = ky
            let atmp = ap[offset - i..]
            for (j in 0..=i) {
                atmp[j] += alpha * (x[ix] * y[jy] + x[jx] * y[iy])
                jx += incX
                jy += incY
            }
            ix += incX
            iy += incY
            offset += i + 2
        }
    }
}