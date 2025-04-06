package cjnum.blas.cjnum

import cjnum.internal.asm.f64.*
import cjnum.blas.*
import std.math.*

extend Float64Implementation <: NFloat64 {
    public func ddot(n: Int64, x: Array<Float64>, incX: Int64, y: Array<Float64>, incY: Int64): Float64 {
        if (incX == 0) {
            throw Exception(zeroIncX)
        }
        if (incY == 0) {
            throw Exception(zeroIncY)
        }
        if (n < 1) {
            if (n == 0) {
                return 0.0
            }
            throw Exception(nLT0)
        }
        if (incX == 1 && incY == 1) {
            if (x.size < n) {
                throw Exception(shortX)
            }
            if (y.size < n) {
                throw Exception(shortY)
            }
            throw Exception("Not implemented") // TODO
        }
        throw Exception("Not implemented")
    }
    
    public func dnrm2(n: Int64, x: Array<Float64>, incX: Int64): Float64 {
        if (incX < 1) {
            if (incX == 0) {
                throw Exception(zeroIncX)
            }
            return 0.0
        }
        if (x.size <= (n - 1) * incX) {
            throw Exception(shortX)
        }
        if (n < 2) {
            match (n) {
                case 1 => return abs(x[0])
                case 0 => return 0.0
                case _ => throw Exception(nLT0)
            }
        }
        if (incX == 1) {
            return l2NormUnitary(x[..n])
        }
        return l2NormInc(x, n, incX)
    }

    public func dasum(n: Int64, x: Array<Float64>, incX: Int64): Float64 {
        var sum: Float64 = 0.0
        if (n < 0) {
            throw Exception(nLT0)
        }
        if (incX < 1) {
            if (incX == 0) {
                throw Exception(zeroIncX)
            }
            return 0.0
        }
        if (x.size <= (n - 1) * incX) {
            throw Exception(shortX)
        }
        if (incX == 1) {
            for (v in x[..n]) {
                sum += abs(v)
            }
            return sum
        }
        for (i in 0..n) {
            sum += abs(x[i * incX])
        }
        return sum
    }

    public func idamax(n: Int64, x: Array<Float64>, incX: Int64): Int64 {
        if (incX < 1) {
            if (incX == 0) {
                throw Exception(zeroIncX)
            }
            return -1
        }
        if (x.size <= (n - 1) * incX) {
            throw Exception(shortX)
        }
        if (n < 2) {
            match (n) {
                case 1 => return 0
                case 0 => return -1
                case _ => throw Exception(nLT0)
            }
        }
        var idx = 0
        var max_ = abs(x[0])
        if (incX == 1) {
            for (i in 0..x[..n].size) {
                let absV = abs(x[i])
                if (absV > max_) {
                    max_ = absV
                    idx = i
                }
            }
            return idx
        }
        var ix = incX
        for (i in 1..n) {
            let absV = abs(x[ix])
            if (absV > max_) {
                max_ = absV
                idx = i
            }
            ix += incX
        }
        return idx
    }

    public func dswap(n: Int64, x: Array<Float64>, incX: Int64, y: Array<Float64>, incY: Int64): Unit {
        if (incX == 0) {
            throw Exception(zeroIncX)
        }
        if (incY == 0) {
            throw Exception(zeroIncY)
        }
        if (n < 1) {
            if (n == 0) {
                return 
            }
            throw Exception(nLT0)
        }
        if ((incX > 0 && x.size <= (n - 1) * incX) || (incX < 0 && x.size <= (1 - n) * incX)) {
            throw Exception(shortX)
        }
        if ((incY > 0 && y.size <= (n - 1) * incY) || (incY < 0 && y.size <= (1 - n) * incY)) {
            throw Exception(shortY)
        }
        if (incX == 1 && incY == 1) {
            let x_ = x[..n]
            for (i in 0..n) {
                let tmp = x_[i]
                x_[i] = y[i]
                y[i] = tmp
            }
            return 
        }
        var ix: Int64 = 0
        var iy: Int64 = 0
        if (incX < 0) {
            ix = (-n + 1) * incX
        }
        if (incY < 0) {
            iy = (-n + 1) * incY
        }
        for (i in 0..n) {
            let tmp = x[ix]
            x[ix] = y[iy]
            y[iy] = tmp
            ix += incX
            iy += incY
        }
    }

    public func dcopy(n: Int64, x: Array<Float64>, incX: Int64, y: Array<Float64>, incY: Int64): Unit {
        if (incX == 0) {
            throw Exception(zeroIncX)
        }
        if (incY == 0) {
            throw Exception(zeroIncY)
        }
        if (n < 1) {
            if (n == 0) {
                return
            }
            throw Exception(nLT0)
        }
        if ((incX > 0 && x.size <= (n - 1) * incX) || (incX < 0 && x.size <= (1 - n) * incX)) {
            throw Exception(shortX)
        }
        if ((incY > 0 && y.size <= (n - 1) * incY) || (incY < 0 && y.size <= (1 - n) * incY)) {
            throw Exception(shortY)
        }
        if (incX == 1 && incY == 1) {
            copy(y[..n], x[..n])
            return 
        }
        var ix: Int64 = 0
        var iy: Int64 = 0
        if (incX < 0) {
            ix = (-n + 1) * incX
        }
        if (incY < 0) {
            iy = (-n + 1) * incY
        }
        for (i in 0..n) {
            y[iy] = x[ix]
            ix += incX
            iy += incY
        }
    }

    public func daxpy(n: Int64, alpha: Float64, x: Array<Float64>, incX: Int64, y: Array<Float64>, incY: Int64): Unit {
        if (incX == 0) {
            throw Exception(zeroIncX)
        }
        if (incY == 0) {
            throw Exception(zeroIncY)
        }
        if (n < 1) {
            if (n == 0) {
                return
            }
            throw Exception(nLT0)
        }
        if ((incX > 0 && x.size <= (n - 1) * incX) || (incX < 0 && x.size <= (1 - n) * incX)) {
            throw Exception(shortX)
        }
        if ((incY > 0 && y.size <= (n - 1) * incY) || (incY < 0 && y.size <= (1 - n) * incY)) {
            throw Exception(shortY)
        }
        if (alpha == 0.0) {
            return 
        }
        if (incX == 1 && incY == 0) {
            axpyUnitary(alpha, x[..n], y[..n])
            return 
        }
        var ix: Int64 = 0
        var iy: Int64 = 0
        if (incX < 0) {
            ix = (-n + 1) * incX
        }
        if (incY < 0) {
            iy = (-n + 1) * incY
        }
        axpyInc(alpha, x, y, n, incX, incY, ix, iy)
    }

    public func drotg(a: Float64, b: Float64): (Float64, Float64, Float64, Float64) {
        const safmin = 0x1p-1022
        const safmax = 1.0 / safmin
        let anorm = abs(a)
        let bnorm = abs(b)
        var c: Float64 = 0.0
        var s: Float64 = 0.0
        var r: Float64 = 0.0
        var z: Float64 = 0.0
        if (bnorm == 0.0) {
            c = 1.0
            s = 0.0
            r = a
            z = 0.0
        } else if (anorm == 0.0) {
            c = 0.0
            s = 1.0
            r = b
            z = 1.0
        } else {
            let maxab = max(anorm, bnorm)
            let scl = min(max(safmin, maxab), safmax)
            var sigma: Float64 = 0.0
            if (anorm > bnorm) {
                sigma = copysign(1.0, a)
            } else {
                sigma = copysign(1.0, b)
            }
            let ascl = a / scl
            let bscl = b / scl
            r = sigma * (scl * sqrt(ascl * ascl + bscl * bscl))
            c = a / r
            s = b / r
            if (anorm > bnorm) {
                z = s
            } else if (c != 0.0) {
                z = 1.0 / c
            } else {
                z = 1.0
            }
        }
        return (c, s ,r ,z)
    }
    
    public func drotmg(d1: Float64, d2: Float64, x1: Float64, y1: Float64): (DrotmParams, Float64, Float64, Float64) {
		const gam    = 4096.0
		const gamsq  = gam * gam
		const rgamsq = 1.0 / gamsq
        let p = DrotmParams()
        if (d1 < 0.0) {
            p.flag = Rescaling
            return (p, 0.0, 0.0, 0.0)
        }

        if (d2 == 0.0 || y1 == 0.0) {
            p.flag = Identity
            return (p, d1, d2, x1)
        }

        var h11: Float64 = 0.0
        var h12: Float64 = 0.0
        var h21: Float64 = 0.0
        var h22: Float64 = 0.0
        var d1_ = d1
        var d2_ = d2
        var x1_ = x1
        if ((d1_ == 0.0 || x1_ == 0.0) && d2_ > 0.0) {
            p.flag = Diagonal
            h12 = 1.0
            h21 = -1.0
            x1_ = y1
            let tmp = d1_
            d1_ = d2_
            d2_ = tmp
        } else {
            let p2 = d2_ * y1
            let p1 = d1_ * x1_
            let q2 = p2 * y1
            let q1 = p1 * x1_
            if (abs(q1) > abs(q2)) {
                p.flag = OffDiagonal
                h11 = 1.0
                h22 = 1.0
                h21 = -y1 / x1_
                h12 = p2 / p1
                let u = 1.0 - (h12 * h21)
                if (u <= 0.0) {
                    p.flag = Rescaling
                    return (p, 0.0, 0.0, 0.0)
                }
                d1_ /= u
                d2_ /= u
                x1_ *= u
            } else {
                if (q2 < 0.0) {
                    p.flag = Rescaling
                    return (p, 0.0, 0.0, 0.0)
                }

                p.flag = Diagonal
                h21 = -1.0
                h12 = 1.0
                h11 = p1 / p2
                h22 = x1_ / y1
                let u = 1.0 + (h11 * h22)
                let new_d1_ = d2_ / u
                let new_d2_ = d1_ / u
                d1_ = new_d1_
                d2_ = new_d2_
                x1_ = y1 * u
            }
        }

        while (d1_ <= rgamsq && d1_ != 0.0) {
            p.flag = Rescaling
            d1_ = (d1_ * gam) * gam
            x1_ /= gam
            h11 /= gam
            h12 /= gam
        }
        while (d1_ > gamsq) {
            p.flag = Rescaling
            d1_ = (d1_ / gam) / gam
            x1_ *= gam
            h11 *= gam
            h12 *= gam
        }

        while (abs(d2_) <= rgamsq && d2_ != 0.0) {
            p.flag = Rescaling
            d2_ = (d2_ * gam) * gam
            h21 /= gam
            h22 /= gam
        }
        while (abs(d2_) > gamsq) {
            p.flag = Rescaling
            d2_ = (d2_ / gam) / gam
            h21 *= gam
            h22 *= gam
        }

        if (p.flag == Diagonal) {
            p.h[0] = h11
            p.h[3] = h22
        } else if (p.flag == OffDiagonal) {
            p.h[1] = h21
            p.h[2] = h12
        } else if (p.flag == Rescaling) {
            p.h = [h11, h21, h12, h22]
        } else {
            throw Exception(badFlag)
        }

        return (p, d1_, d2_, x1_)
    }

    public func drot(n: Int64, x: Array<Float64>, incX: Int64, y: Array<Float64>, incY: Int64, c: Float64, s: Float64): Unit {
        if (incX == 0) {
            throw Exception(zeroIncX)
        }
        if (incY == 0) {
            throw Exception(zeroIncY)
        }
        if (n < 1) {
            if (n == 0) {
                return 
            }
            throw Exception(nLT0)
        }
        if ((incX > 0 && x.size <= (n - 1) * incX) || (incX < 0 && x.size <= (1 - n) * incX)) {
            throw Exception(shortX)
        }
        if ((incY > 0 && y.size <= (n - 1) * incY) || (incY < 0 && y.size <= (1 - n) * incY)) {
            throw Exception(shortY)
        }
        if (incX == 1 && incY == 1) {
            let x_ = x[..n]
            for (i in 0..n) { 
                let vy = y[i]
                let new_vx = c * x_[i] + s * vy
                let new_vy = c * vy - s * x_[i]
                x_[i] = new_vx
                y[i] = new_vy
            }
            return
        }
        var ix = 0
        var iy = 0
        if (incX < 0) {
            ix = (-n + 1) * incX
        }
        if (incY < 0) {
            iy = (-n + 1) * incY
        }
        for (i in 0..n) {
            let vx = x[ix]
            let vy = y[iy]
            let new_vx = c * vx + s * vy
            let new_vy = c * vy - s * vx
            x[ix] = new_vx
            y[iy] = new_vy
            ix += incX
            iy += incY
        }
    }

    public func drotm(n: Int64, x: Array<Float64>, incX: Int64, y: Array<Float64>, incY: Int64, p: DrotmParams): Unit {
        if (incX == 0) {
            throw Exception(zeroIncX)
        }
        if (incY == 0) {
            throw Exception(zeroIncY)
        }
        if (n <= 0) {
            if (n == 0) {
                return 
            }
            throw Exception(nLT0)
        }
        if ((incX > 0 && x.size <= (n - 1) * incX) || (incX < 0 && x.size <= (1 - n) * incX)) {
            throw Exception(shortX)
        }
        if ((incY > 0 && y.size <= (n - 1) * incY) || (incY < 0 && y.size <= (1 - n) * incY)) {
            throw Exception(shortY)
        }

        if (p.flag == Identity) {
            return 
        }

        if (p.flag == Rescaling) {
            let h11 = p.h[0]
            let h12 = p.h[2]
            let h21 = p.h[1]
            let h22 = p.h[3]
            if (incX == 1 && incY == 1) {
                for (i in 0..n) {
                    let vy = y[i]
                    let new_vx = h11 * x[i] + h12 * vy
                    let new_vy = h21 * x[i] + h22 * vy
                    x[i] = new_vx
                    y[i] = new_vy
                }
                return 
            }
            var ix = 0
            var iy = 0
            if (incX < 0) {
                ix = (-n + 1) * incX
            }
            if (incY < 0) {
                iy = (-n + 1) * incY
            }
            for (i in 0..n) {
                let vx = x[ix]
                let vy = y[iy]
                let new_vx = h11 * vx + h12 * vy
                let new_vy = h21 * vx + h22 * vy
                x[ix] = new_vx
                y[iy] = new_vy
                ix += incX
                iy += incY
            }
        }
        else if (p.flag == OffDiagonal) {
            let h12 = p.h[2]
            let h21 = p.h[1]
            if (incX == 1 && incY == 1) {
                let x_ = x[..n]
                for (i in 0..x_.size) {
                    let vy = y[i]
                    let new_vx = x_[i] + h12 * vy
                    let new_vy = h21 * x_[i] + vy
                    x_[i] = new_vx
                    y[i] = new_vy
                }
                return 
            }
            var ix = 0
            var iy = 0
            if (incX < 0) {
                ix = (-n + 1) * incX
            }
            if (incY < 0) {
                iy = (-n + 1) * incY
            }
            for (i in 0..n) {
                let vx = x[ix]
                let vy = y[iy]
                let new_vx = vx + h12 * vy
                let new_vy = h21 * vx + vy
                x[ix] = new_vx
                y[iy] = new_vy
                ix += incX
                iy += incY
            }
        }
        else if (p.flag == Diagonal) {
            let h11 = p.h[0]
            let h22 = p.h[3]
            if (incX == 1 && incY == 1) {
                for (i in 0..n) {
                    let vy = y[i]
                    let new_vx = h11 * x[i] + vy
                    let new_vy = -x[i] + h22 * vy
                    x[i] = new_vx
                    y[i] = new_vy
                }
                return 
            }
            var ix = 0
            var iy = 0
            if (incX < 0) {
                ix = (-n + 1) * incX
            }
            if (incY < 0) {
                iy = (-n + 1) * incY
            }
            for (i in 0..n) {
                let vx = x[ix]
                let vy = y[iy]
                let new_vx = h11 * vx + vy
                let new_vy = -vx + h22 * vy
                x[ix] = new_vx
                y[iy] = new_vy
                ix += incX
                iy += incY
            }
        }

        return 
    }

    public func dscal(n: Int64, alpha: Float64, x: Array<Float64>, incX: Int64): Unit {
        if (incX < 1) {
            if (incX == 0) {
                throw Exception(zeroIncX)
            }
            return 
        }
        if (n < 1) {
            if (n == 0) {
                return 
            }
            throw Exception(nLT0)
        }
        if (alpha == 0.0) {
            if (incX == 1) {
                let x_ = x[..n]
                for (i in 0..x_.size) {
                    x_[i] = 0.0
                }
                return 
            }
            var ix = 0
            while (ix < n * incX) {
                x[ix] = 0.0
                ix += incX
            }
            return 
        }
        if (incX == 1) {
            scalUnitary(alpha, x[..n])
            return 
        }
        scalInc(alpha, x, n, incX)
        return 
    }
}

func copy<T>(dst: Array<T>, src: Array<T>): Int64 {
    let minSize = min(src.size, dst.size)
    dst[0..minSize] = src[0..minSize]
    return minSize
}

func copysign(f: Float64, sign: Float64): Float64 {
    const signBit: UInt64 = 1 << 63
    return Float64.fromBits((f.toBits() & (!signBit)) | sign.toBits() & signBit)
}
