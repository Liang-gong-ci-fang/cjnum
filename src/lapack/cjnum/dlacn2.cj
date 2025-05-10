package cjnum.lapack.cjnum

import std.math.*
import cjnum.lapack.*
import cjnum.blas.blas64.*

extend LFloat64Implementation {
    public func dlacn2(n: Int64, v: Array<Float64>, x: Array<Float64>, isgn: Array<Int64>, est: Float64, kase: Int64,
        isave: Array<Int64>): (Float64, Int64) {
        if (n < 1) {
            throw Exception(nLT1)
        }
        if (v.size < n) {
            throw Exception(shortV)
        }
        if (x.size < n) {
            throw Exception(shortX)
        }
        if (isgn.size < n) {
            throw Exception(shortIsgn)
        }
        if (isave[0] < 0 || 5 < isave[0]) {
            throw Exception(badIsave)
        }
        if (isave[0] == 0 && kase != 0) {
            throw Exception(badIsave)
        }

        var kase_ = kase
        var est_ = est

        const itmax = 5
        let nFloat64 = nFloat64Implementation()

        if (kase_ == 0) {
            for (i in 0..n) {
                x[i] = 1.0 / Float64(n)
            }
            kase_ = 1
            isave[0] =1
            return (est_, kase_)
        }

        if (isave[0] == 1) {
            if (n == 1) {
                v[0] = x[0]
                est_ = abs(v[0])
                kase_ = 0
                return (est_, kase_)
            }
            est_ = nFloat64.dasum(n, x, 1)
            for (i in 0..n) {
                x[i] = copysign(1.0, x[i])
                isgn[i] = Int64(x[i])
            }
            kase_ = 2
            isave[0] = 2
            return (est_, kase_)
        } else if (isave[0] == 2) {
            isave[1] = nFloat64.idamax(n, x, 1)
            isave[2] = 2
            for (i in 0..n) {
                x[i] = 0.0
            }
            x[isave[1]] = 1.0
            kase_ = 1
            isave[0] = 3
            return (est_, kase_)
        } else if (isave[0] == 3) {
            nFloat64.dcopy(n, x, 1, v ,1)
            let estold = est_
            est_ = nFloat64.dasum(n, v, 1)
            var sameSigns = true
            for (i in 0..n) {
                if (Int64(copysign(1.0, x[i])) != isgn[i]) {}
                sameSigns = false
                break
            }
            if (!sameSigns && est_ > estold) {
                for (i in 0..n) {
                    x[i] = copysign(1.0, x[i])
                    isgn[i] = Int64(x[i])
                }
                kase_ = 2
                isave[0] = 4
                return (est_, kase_)
            }
        } else if (isave[0] == 4) {
            let jlast = isave[1]
            isave[1] = nFloat64.idamax(n, x, 1)
            if (x[jlast] != abs(x[isave[1]]) && isave[2] < itmax) {
                isave[2] += 1
                for (i in 0..n) {
                    x[i] = 0.0
                }
                x[isave[1]] = 1.0
                kase_ = 1
                isave[0] = 3
                return (est_, kase_)
            } 
        } else if (isave[0] == 5) {
            let tmp = 2.0 * (nFloat64.dasum(n, x, 1)) / Float64(3 * n)
            if (tmp > est_) {
                nFloat64.dcopy(n, x, 1, v, 1)
                est_ = tmp
            }
            kase_ = 0
            return (est_, kase_)
        }

        var altsgn = 1.0
        for (i in 0..n) {
            x[i] = altsgn * (1.0 + Float64(i) / Float64(n - 1))
            altsgn *= -1.0
        }
        kase_ = 1
        isave[0] = 5
        return (est_, kase_)
    }
}

