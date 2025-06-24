package cjnum.lapack.cjnum

import cjnum.blas.blas64.*
import std.math.*

extend LFloat64Implementation {
    public func dlapmr(forward: Bool, m: Int64, n: Int64, x: Array<Float64>, ldx: Int64, k: Array<Int64>): Unit {
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
        if (k.size != m) {
            throw Exception(badLenK)
        }

        if (m == 1) {
            return
        }

        let nFloat64 = nFloat64Implementation()

        for (i in 0..k.size) {
            let ki = k[i]
            k[i] = -(ki + 1)
        }
        if (forward) {
            for (i in 0..k.size) {
                let ki = k[i]
                if (ki >= 0) {
                    continue
                }
                var j = i
                k[j] = -k[j] - 1
                var in_ = k[j]
                while (true) {
                    if (k[in_] >= 0) {
                        break
                    }
                    nFloat64.dswap(n, x[j * ldx..], 1, x[in_ * ldx..], 1)
                    k[in_] = -k[in_] - 1
                    j = in_
                    in_ = k[in_]
                }
            }
        } else {
            for (i in 0..k.size) {
                let ki = k[i]
                if (ki >= 0) {
                    continue
                }
                k[i] = -ki - 1
                var j = k[i]
                while (true) {
                    if (j == i) {
                        break
                    }
                    nFloat64.dswap(n, x[i * ldx..], 1, x[j * ldx..], 1)
                    k[j] = -k[j] - 1
                    j = k[j]
                }
            }
        }
    }
}
