package cjnum.lapack.cjnum

import std.math.*
import std.math.numeric.*
import cjnum.blas.blas64.*

extend LFloat64Implementation {
    public func dlarfg(n: Int64, alpha: Float64, x: Array<Float64>, incX: Int64): (Float64, Float64) {
        if (n < 0) {
            throw Exception(nLT0)
        }
        if (incX <= 0) {
            throw Exception(badIncX)
        }

        if (n <= 1) {
            return (alpha, 0.0)
        }

        if (x.size < 1 + (n - 2) * abs(incX)) {
            throw Exception(shortX)
        }

        let nFloat64 = nFloat64Implementation()

        var xnorm = nFloat64.dnrm2(n - 1, x, incX)
        var alpha_ = alpha
        if (xnorm == 0.0) {
            return (alpha_, 0.0)
        }
        var beta = -copysign(dlapy2(alpha_, xnorm), alpha_)
        let safmin = dlamchS / dlamchE
        var knt = 0
        if (abs(beta) < safmin) {
            let rsafmn = 1.0 / safmin
            while (true) {
                knt++
                nFloat64.dscal(n - 1, rsafmn, x, incX)
                beta *= rsafmn
                alpha_ *= rsafmn
                if (abs(beta) >= safmin) {
                    break
                }
            }
            xnorm = nFloat64.dnrm2(n - 1, x, incX)
            beta = -copysign(dlapy2(alpha_, xnorm), alpha_)
        }
        let tau = (beta - alpha_) / beta
        nFloat64.dscal(n - 1, 1.0 / (alpha_ - beta), x, incX)
        for (j in 0..knt) {
            beta *= safmin
        }
        return (beta, tau)
    }
}