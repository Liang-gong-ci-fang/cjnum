package cjnum.lapack.cjnum

import std.math.*
import std.math.numeric.*

extend LFloat64Implementation {
    public func dlaqr1(n: Int64, h: Array<Float64>, ldh: Int64, sr1: Float64, si1: Float64, sr2: Float64, si2: Float64,
        v: Array<Float64>): Unit {
        if (n != 2 && n != 3) {
            throw Exception("lapack: n must be 2 or 3")
        }
        if (ldh < n) {
            throw Exception(badLdH)
        }
        if (h.size < (n - 1) * ldh + n) {
            throw Exception(shortH)
        }
        if (!((sr1 == sr2 && si1 == -si2) || (si1 == 0.0 && si2 == 0.0))) {
            throw Exception(badShifts)
        }
        if (v.size != n) {
            throw Exception(shortV)
        }

        if (n == 2) {
            let s = abs(h[0] - sr2) + abs(si2) + abs(h[ldh])
            if (s == 0.0) {
                v[0] = 0.0
                v[1] = 0.0
            } else {
                let h21s = h[ldh] / s
                v[0] = h21s * h[1] + (h[0] - sr1) * ((h[0] - sr2) / s) - si1 * (si2 / s)
                v[1] = h21s * (h[0] + h[ldh + 1] - sr1 - sr2)
            }
            return
        }

        let s = abs(h[0] - sr2) + abs(si2) + abs(h[ldh]) + abs(h[2 * ldh])
        if (s == 0.0) {
            v[0] = 0.0
            v[1] = 0.0
            v[2] = 0.0
        } else {
            let h21s = h[ldh] / s
            let h31s = h[2 * ldh] / s
            v[0] = (h[0] - sr1) * ((h[0] - sr2) / s) - si1 * (si2 / s) + h[1] * h21s + h[2] * h31s
            v[1] = h21s * (h[0] + h[ldh + 1] - sr1 - sr2) + h[ldh + 2] * h31s
            v[2] = h31s * (h[0] + h[2 * ldh + 2] - sr1 - sr2) + h21s * h[2 * ldh + 1]
        }
    }
}
