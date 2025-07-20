package cjnum.lapack.cjnum

import cjnum.blas.blas64.*
import std.math.*

extend LFloat64Implementation {
    public func dlapmt(forward: Bool, m: Int64, n: Int64, x: Array<Float64>, ldx: Int64, k: Array<Int64>): Unit {
        if (m < 0) {
            throw Exception(mLT0)
        }
        if (n < 0) {
            throw Exception(nLT0)
        }
        if (ldx < max(1, n)) {
            throw Exception(badLdX)
        }

        if (m == 0 || n == 0) {
            return
        }

        if (x.size < (m - 1) * ldx + n) {
            throw Exception(shortX)
        }
        if (k.size != n) {
            throw Exception(badLenK)
        }

        if (n == 1) {
            return
        }

        for (i in 0..k.size) {
            var v = k[i]
            v++
            k[i] = -v
        }

        let nFloat64 = nFloat64Implementation()

        if (forward) {
            var j = 0
            while (j < k.size) {
                let v = k[j]
                if (v >= 0) {
                    j++
                    continue
                }
                k[j] = -v
                var i = -v - 1
                while (k[i] < 0) {
                    nFloat64.dswap(m, x[j..], ldx, x[i..], ldx)

                    k[i] = -k[i]
                    j = i
                    i = k[i] - 1
                }
                j++
            }
        } else {
            for (i in 0..k.size) {
                let v = k[i]
                if (v >= 0) {
                    continue
                }
                k[i] = -v
                var j = -v - 1
                while (j != i) {
                    nFloat64.dswap(m, x[j..], ldx, x[i..], ldx)

                    k[j] = -k[j]
                    j = k[j] - 1
                }
            }
        }

        for (i in 0..k.size) {
            k[i]--
        }
    }
}
