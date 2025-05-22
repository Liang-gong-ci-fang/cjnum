package cjnum.lapack.cjnum

import std.math.*
import cjnum.lapack.*

extend LFloat64Implementation {
    public func dlanst(norm: MatrixNorm, n: Int64, d: Array<Float64>, e: Array<Float64>): Float64 {
        if (norm != MaxAbs && norm != MaxRowSum && norm != MaxColumnSum && norm != Frobenius) {
            throw Exception(badNorm)
        }
        if (n < 0) {
            throw Exception(nLT0)
        }
        if (n == 0) {
            return 0.0
        }
        if (d.size < n) {
            throw Exception(shortD)
        }
        if (e.size < n - 1) {
            throw Exception(shortE)
        }

        match {
            case norm == MaxAbs =>
                var anorm = abs(d[n - 1])
                for (i in 0..n - 1) {
                    var sum = abs(d[i])
                    if (anorm < sum || sum.isNaN()) {
                        anorm = sum
                    }
                    sum = abs(e[i])
                    if (anorm < sum || sum.isNaN()) {
                        anorm = sum
                    }
                }
                return anorm
            case norm == MaxColumnSum || norm == MaxRowSum =>
                if (n == 1) {
                    return abs(d[0])
                }
                var anorm = abs(d[0]) + abs(e[0])
                var sum = abs(e[n - 2]) + abs(d[n - 1])
                if (anorm < sum || sum.isNaN()) {
                    anorm = sum
                }
                for (i in 1..n - 1) {
                    sum = abs(d[i]) + abs(e[i]) + abs(e[i-1])
                    if (anorm < sum || sum.isNaN()) {
                        anorm = sum
                    }
                }
                return anorm
            case norm == Frobenius =>
                var scale = 0.0
                var sum = 1.0
                if (n > 1) {
                    (scale, sum) = dlassq(n - 1, e, 1, scale, sum)
                    sum = 2.0 * sum
                }
                (scale, sum) = dlassq(n, d, 1, scale, sum)
                return scale * sqrt(sum)
            case _ => throw Exception(badNorm)
        }
    }
}