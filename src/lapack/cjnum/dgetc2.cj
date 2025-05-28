package cjnum.lapack.cjnum

import cjnum.blas.blas64.*
import cjnum.lapack.*
import std.math.*

extend LFloat64Implementation {
    public func dgetc2(n: Int, a: Array<Float64>, lda: Int, ipiv: Array<Int>, jpiv: Array<Int>): Int {
        if (n < 0) {
            throw Exception(nLT0)
        } else if (lda < max(1, n)) {
            throw Exception(badLdA)
        }

        var k = -1

        if (n == 0) {
            return k
        }

        if (a.size < (n - 1) * lda + n) {
            throw Exception(shortA)
        } else if (ipiv.size != n) {
            throw Exception(badLenIpiv)
        } else if (jpiv.size != n) {
            throw Exception(badLenJpvt)
        }

        const eps = dlamchP
        const smlnum = dlamchS / eps

        if (n == 1) {
            (ipiv[0], jpiv[0]) = (0, 0)
            if (abs(a[0]) < smlnum) {
                k = 0
                a[0] = smlnum
            }
            return k
        }

        var smin: Float64 = 0.0
        var (ipv, jpv) = (0, 0)
        let nFloat64 = nFloat64Implementation()

        for (i in 0..n - 1) {
            var xmax: Float64 = 0.0
            for (ip in i..n) {
                for (jp in i..n) {
                    if (abs(a[ip * lda + jp]) >= xmax) {
                        xmax = abs(a[ip * lda + jp])
                        ipv = ip
                        jpv = jp
                    }
                }
            }
            if (i == 0) {
                smin = max(eps * xmax, smlnum)
            }

            if (ipv != i) {
                nFloat64.dswap(n, a[ipv * lda..], 1, a[i * lda..], 1)
            }
            ipiv[i] = ipv

            if (jpv != i) {
                nFloat64.dswap(n, a[jpv..], lda, a[i..], lda)
            }
            jpiv[i] = jpv

            if (abs(a[i * lda + i]) < smin) {
                k = i
                a[i * lda + i] = smin
            }

            for (j in i + 1..n) {
                a[j * lda + i] /= a[i * lda + i]
            }

            nFloat64.dger(n - i - 1, n - i - 1, -1.0, a[(i + 1) * lda + i..], lda, a[i * lda + i + 1..], 1,
                a[(i + 1) * lda + i + 1..], lda)
        }

        if (abs(a[(n - 1) * lda + (n - 1)]) < smin) {
            k = n - 1
            a[(n - 1) * lda + (n - 1)] = smin
        }

        ipiv[n - 1] = n - 1
        jpiv[n - 1] = n - 1

        return k
    }
}
