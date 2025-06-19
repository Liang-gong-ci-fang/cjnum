package cjnum.lapack.cjnum

import cjnum.lapack.*
import std.math.*

extend LFloat64Implementation {
    public func dlangt(norm: MatrixNorm, n: Int64, dl: Array<Float64>, d: Array<Float64>, du: Array<Float64>): Float64 {
        var dl_: Array<Float64> = dl
        var d_: Array<Float64> = d
        var du_: Array<Float64> = du

        if (norm != MaxAbs && norm != MaxRowSum && norm != MaxColumnSum && norm != Frobenius) {
            throw Exception(badNorm)
        }
        if (n < 0) {
            throw Exception(nLT0)
        }

        if (n == 0) {
            return 0.0
        }

        if (dl_.size < n - 1) {
            throw Exception(shortDL)
        }
        if (d_.size < n) {
            throw Exception(shortD)
        }
        if (du_.size < n - 1) {
            throw Exception(shortDU)
        }

        dl_ = dl_[0..n - 1]
        d_ = d_[0..n]
        du_ = du_[0..n - 1]

        var anorm = 0.0
        if (norm == MaxAbs) {
            for (diag in [dl_, d_, du_]) {
                for (i in 0..diag.size) {
                    var di = diag[i]
                    if (di.isNaN()) {
                        return di
                    }
                    di = abs(di)
                    if (di > anorm) {
                        anorm = di
                    }
                }
            }
        } else if (norm == MaxColumnSum) {
            if (n == 1) {
                return abs(d_[0])
            }
            anorm = abs(d_[0]) + abs(dl_[0])
            if (anorm.isNaN()) {
                return anorm
            }
            var tmp = abs(du_[n - 2]) + abs(d_[n - 1])
            if (tmp.isNaN()) {
                return tmp
            }
            if (tmp > anorm) {
                anorm = tmp
            }
            for (i in 1..n - 1) {
                tmp = abs(du_[i - 1]) + abs(d_[i]) + abs(dl_[i])
                if (tmp.isNaN()) {
                    return tmp
                }
                if (tmp > anorm) {
                    anorm = tmp
                }
            }
        } else if (norm == MaxRowSum) {
            if (n == 1) {
                return abs(d_[0])
            }
            anorm = abs(d_[0]) + abs(du_[0])
            if (anorm.isNaN()) {
                return anorm
            }
            var tmp = abs(dl_[n - 2]) + abs(d_[n - 1])
            if (tmp.isNaN()) {
                return tmp
            }
            if (tmp > anorm) {
                anorm = tmp
            }
            for (i in 1..n - 1) {
                tmp = abs(dl_[i - 1]) + abs(d_[i]) + abs(du_[i])
                if (tmp.isNaN()) {
                    return tmp
                }
                if (tmp > anorm) {
                    anorm = tmp
                }
            }
        } else if (norm == Frobenius) {
            var scale = 0.0
            var ssq = 1.0
            (scale, ssq) = dlassq(n, d_, 1, scale, ssq)
            if (n > 1) {
                (scale, ssq) = dlassq(n - 1, dl_, 1, scale, ssq)
                (scale, ssq) = dlassq(n - 1, du_, 1, scale, ssq)
            }
            anorm = scale * sqrt(ssq)
        }
        return anorm
    }
}
