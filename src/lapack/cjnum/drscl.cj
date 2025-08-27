package cjnum.lapack.cjnum

import std.math.*
import cjnum.blas.blas64.*

extend LFloat64Implementation {
    public func drscl(n: Int64, a: Float64, x: Array<Float64>, incX: Int64): Unit {
        if (n < 0) {
            throw Exception(nLT0)
        }
        if (incX <= 0) {
            throw Exception(badIncX)
        }

        if (n == 0) {
            return 
        }

        if (x.size < 1 + (n - 1) * incX) {
            throw Exception(shortX)
        }

        let nFloat64 = nFloat64Implementation()

        var cden = a
        var cnum = 1.0
        let smlnum = dlamchS
        let bignum = 1.0 / smlnum
        while (true) {
            let cden1 = cden * smlnum
            let cnum1 = cnum  / bignum
            var mul = 0.0
            var done = false
            if (cnum != 0.0 && abs(cden1) > abs(cnum)) {
                mul = smlnum
                done = false
                cden = cden1
            } else if (abs(cnum1) > abs(cden)) {
                mul = bignum
                done = false
                cnum = cnum1
            } else {
                mul = cnum / cden
                done = true
            }
            nFloat64.dscal(n, mul, x, incX)
            if (done) {
                break
            }
        }
    }
}