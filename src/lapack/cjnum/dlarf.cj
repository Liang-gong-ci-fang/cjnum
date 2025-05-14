package cjnum.lapack.cjnum

import std.math.*
import cjnum.blas.*
import cjnum.blas.blas64.*

extend LFloat64Implementation {
    public func dlarf(side: Side, m: Int64, n: Int64, v: Array<Float64>, incv: Int64, tau: Float64, c: Array<Float64>,
        ldc: Int64, work: Array<Float64>): Unit {
        if (side != Left && side != Right) {
            throw Exception(badSide)
        }
        if (m < 0) {
            throw Exception(mLT0)
        }
        if (n < 0) {
            throw Exception(nLT0)
        }
        if (incv == 0) {
            throw Exception(zeroIncV)
        }
        if (ldc < max(1, n)) {
            throw Exception(badLdC)
        }

        if (m == 0 || n == 0) {
            return
        }

        let applyleft = side == Left
        let lenV = if (applyleft) {
            m
        } else {
            n
        }

        if (v.size < 1 + (lenV - 1) * abs(incv)) {
            throw Exception(shortV)
        }
        if (c.size < (m - 1) * ldc + n) {
            throw Exception(shortC)
        }
        if ((applyleft && work.size < n) || (!applyleft && work.size < m)) {
            throw Exception(shortWork)
        }

        var lastv = -1
        var lastc = -1
        if (tau != 0.0) {
            if (applyleft) {
                lastv = m - 1
            } else {
                lastv = n - 1
            }
            var i = 0
            if (incv > 0) {
                i = lastv * incv
            }
            while (lastv >= 0 && v[i] == 0.0) {
                lastv--
                i -= incv
            }
            if (applyleft) {
                lastc = iladlc(lastv + 1, n, c, ldc)
            } else {
                lastc = iladlr(m, lastv + 1, c, ldc)
            }
        }
        if (lastv == -1 || lastc == -1) {
            return
        }
        let nFloat64 = nFloat64Implementation()
        if (applyleft) {
            nFloat64.dgemv(Trans, lastv + 1, lastc + 1, 1.0, c, ldc, v, incv, 0.0, work, 1)
            nFloat64.dger(lastv + 1, lastc + 1, -tau, v, incv, work, 1, c, ldc)
        } else {
            nFloat64.dgemv(NoTrans, lastc + 1, lastv + 1, 1.0, c, ldc, v, incv, 0.0, work, 1)
            nFloat64.dger(lastc + 1, lastv + 1, -tau, work, 1, v, incv, c, ldc)
        }
    }
}
