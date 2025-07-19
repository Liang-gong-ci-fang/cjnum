package cjnum.lapack.cjnum

import cjnum.blas.blas64.*
import cjnum.blas.*
import cjnum.lapack.General as lGeneral
import cjnum.lapack.*
import std.math.*

extend LFloat64Implementation {
    public func dlasq1(n: Int64, d: Array<Float64>, e: Array<Float64>, work: Array<Float64>): Int64 {
        var info = 0
        if (n < 0) {
            throw Exception(nLT0)
        }

        if (n == 0) {
            return info
        }

        if (d.size < n) {
            throw Exception(shortD)
        } else if (e.size < n - 1) {
            throw Exception(shortE)
        } else if (work.size < 4 * n) {
            throw Exception(shortWork)
        }

        if (n == 1) {
            d[0] = abs(d[0])
            return info
        }

        if (n == 2) {
            (d[1], d[0]) = dlas2(d[0], e[0], d[1])
            return info
        }

        var sigmx: Float64 = 0.0
        for (i in 0..n - 1) {
            d[i] = abs(d[i])
            sigmx = max(sigmx, abs(e[i]))
        }
        d[n - 1] = abs(d[n - 1])
        if (sigmx == 0.0) {
            dlasrt(SortDecreasing, n, d)
            return info
        }

        for (i in 0..n) {
            sigmx = max(sigmx, d[i])
        }

        let eps = dlamchP
        let safmin = dlamchS
        let scale = sqrt(eps / safmin)
        let nFloat64 = nFloat64Implementation()
        nFloat64.dcopy(n, d, 1, work, 2)
        nFloat64.dcopy(n - 1, e, 1, work[1..], 2)
        dlascl(lGeneral, 0, 0, sigmx, scale, 2 * n - 1, 1, work, 1)

        for (i in 0..2 * n - 1) {
            work[i] *= work[i]
        }
        work[2 * n - 1] = 0.0

        info = dlasq2(n, work)
        if (info == 0) {
            for (i in 0..n) {
                d[i] = sqrt(work[i])
            }
            dlascl(lGeneral, 0, 0, scale, sigmx, n, 1, d, 1)
        } else if (info == 2) {
            for (i in 0..n) {
                d[i] = sqrt(work[2 * i])
                e[i] = sqrt(work[2 * i + 1])
            }
            dlascl(lGeneral, 0, 0, scale, sigmx, n, 1, d, 1)
            dlascl(lGeneral, 0, 0, scale, sigmx, n, 1, e, 1)
        }
        return info
    }
}
