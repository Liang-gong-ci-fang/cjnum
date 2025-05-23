package cjnum.lapack.cjnum

import std.math.*
import cjnum.lapack.*
import cjnum.blas.*
import cjnum.blas.blas64.*

extend LFloat64Implementation {
    public func dgetf2(m: Int64, n: Int64, a: Array<Float64>, lda: Int64, ipiv: Array<Int64>): Bool {
        let mn = min(m, n)
        if (m < 0) {
            throw Exception(mLT0)
        }
        if (n < 0) {
            throw Exception(nLT0)
        }
        if (lda < max(1, n)) {
            throw Exception(badLdA)
        }

        if (mn == 0) {
            return true
        }

        if (a.size < (m - 1) * lda + n) {
            throw Exception(shortA)
        } 
        if (ipiv.size != mn) {
            throw Exception(badLenIpiv)
        }

        let nFloat64 = nFloat64Implementation()
        let sfmin = dlamchS
        var ok = true
        for (j in 0..mn) {
            let jp = j + nFloat64.idamax(m - j, a[j * lda + j..], lda)
            ipiv[j] = jp
            if (a[jp * lda + j] == 0.0) {
                ok = false
            } else {
                if (jp != j) {
                    nFloat64.dswap(n, a[j * lda..], 1, a[jp * lda..], 1)
                }
                if (j < m - 1) {
                    let aj = a[j * lda + j]
                    if (abs(aj) >= sfmin) {
                        nFloat64.dscal(m - j - 1, 1.0 / aj, a[(j + 1) * lda + j..], lda)
                    } else {
                        for (_ in 0..m - j - 1) {
                            a[(j + 1) * lda + j] = a[(j + 1) * lda + j] / a[lda * j + j]
                        }
                    }
                }
            }
            if (j < mn - 1) {
                nFloat64.dger(m - j - 1, n - j - 1, -1.0, a[(j + 1) * lda + j..], lda, a[j * lda + j + 1..], 1, a[(j + 1) * lda + j + 1..], lda)
            }
        }
        return ok
    }
}