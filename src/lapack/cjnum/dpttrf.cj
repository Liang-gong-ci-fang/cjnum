package cjnum.lapack.cjnum

import std.math.*
import std.math.numeric.*

extend LFloat64Implementation {
    public func dpttrf(n: Int64, d: Array<Float64>, e: Array<Float64>): Bool {
        if (n < 0) {
            throw Exception(nLT0)
        }
        if (n == 0) {
            return true
        }

        if (d.size < n) {
            throw Exception(shortD)
        }
        if (e.size < n - 1) {
            throw Exception(shortE)
        }

        let i4 = (n - 1) % 4
        for (i in 0..i4) {
            if (d[i] <= 0.0) {
                return false
            }
            let ei = e[i]
            e[i] /= d[i]
            d[i + 1] -= e[i] * ei
        }

        var i = i4
        while (i < n - 4) {
            if (d[i] <= 0.0) {
                return false
            }

            var ei = e[i]
            e[i] /= d[i]
            d[i + 1] -= e[i] * ei
            if (d[i + 1] <= 0.0) {
                return false
            }

            ei = e[i+1]
            e[i + 1] /= d[i + 1]
            d[i + 2] -= e[i + 1] * ei
            if (d[i + 2] <= 0.0) {
                return false
            }

            // Solve for e[i+2] and d[i+3].
            ei = e[i + 2]
            e[i + 2] /= d[i + 2]
            d[i + 3] -= e[i + 2] * ei
            if (d[i + 3] <= 0.0) {
                return false
            }

            // Solve for e[i+3] and d[i+4].
            ei = e[i + 3]
            e[i + 3] /= d[i + 3]
            d[i + 4] -= e[i + 3] * ei
            i += 4
        }
        return d[n - 1] > 0.0
    }
}