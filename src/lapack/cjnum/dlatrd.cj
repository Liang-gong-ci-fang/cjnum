package cjnum.lapack.cjnum

import cjnum.blas.blas64.*
import cjnum.blas.*

extend LFloat64Implementation {
    public func dlatrd(uplo: Uplo, n: Int64, nb: Int64, a: Array<Float64>, lda: Int64, e: Array<Float64>, tau: Array<Float64>,
        w: Array<Float64>, ldw: Int64): Unit {
        if (uplo != Upper && uplo != Lower) {
            throw Exception(badUplo)
        } else if (n < 0) {
            throw Exception(nLT0)
        } else if (nb < 0) {
            throw Exception(nbLT0)
        } else if (nb > n) {
            throw Exception(nbGTN)
        } else if (lda < max(1, n)) {
            throw Exception(badLdA)
        } else if (ldw < max(1, nb)) {
            throw Exception(badLdW)
        }

        if (n == 0) {
            return
        }

        if (a.size < (n - 1) * lda + n) {
            throw Exception(shortA)
        } else if (w.size < (n - 1) * ldw + nb) {
            throw Exception(shortW)
        } else if (e.size < n - 1) {
            throw Exception(shortE)
        } else if (tau.size < n - 1) {
            throw Exception(shortTau)
        }

        let nFloat64 = nFloat64Implementation()

        if (uplo == Upper) {
            for (i in n - 1..=n - nb : -1) {
                let iw = i - n + nb
                if (i < n - 1) {
                    nFloat64.dgemv(NoTrans, i + 1, n - i - 1, -1.0, a[i + 1..], lda, w[i * ldw + iw + 1..], 1, 1.0,
                        a[i..], lda)
                    nFloat64.dgemv(NoTrans, i + 1, n - i - 1, -1.0, w[iw + 1..], ldw, a[i * lda + i + 1..], 1, 1.0,
                        a[i..], lda)
                }
                if (i > 0) {
                    (e[i - 1], tau[i - 1]) = dlarfg(i, a[(i - 1) * lda + i], a[i..], lda)
                    a[(i - 1) * lda + i] = 1.0
                    nFloat64.dsymv(Upper, i, 1.0, a, lda, a[i..], lda, 0.0, w[iw..], ldw)
                    if (i < n - 1) {
                        nFloat64.dgemv(Trans, i, n - i - 1, 1.0, w[iw + 1..], ldw, a[i..], lda, 0.0,
                            w[(i + 1) * ldw + iw..], ldw)
                        nFloat64.dgemv(NoTrans, i, n - i - 1, -1.0, a[i + 1..], lda, w[(i + 1) * ldw + iw..], ldw, 1.0,
                            w[iw..], ldw)
                        nFloat64.dgemv(Trans, i, n - i - 1, 1.0, a[i + 1..], lda, a[i..], lda, 0.0, w[(i + 1) * ldw + iw..],
                            ldw)
                        nFloat64.dgemv(NoTrans, i, n - i - 1, -1.0, w[iw + 1..], ldw, w[(i + 1) * ldw + iw..], ldw, 1.0,
                            w[iw..], ldw)
                    }
                    nFloat64.dscal(i, tau[i - 1], w[iw..], ldw)
                    let alpha = -0.5 * tau[i - 1] * nFloat64.ddot(i, w[iw..], ldw, a[i..], lda)
                    nFloat64.daxpy(i, alpha, a[i..], lda, w[iw..], ldw)
                }
            }
        } else {
            for (i in 0..nb) {
                nFloat64.dgemv(NoTrans, n - i, i, -1.0, a[i * lda..], lda, w[i * ldw..], 1, 1.0, a[i * lda + i..], lda)
                nFloat64.dgemv(NoTrans, n - i, i, -1.0, w[i * ldw..], ldw, a[i * lda..], 1, 1.0, a[i * lda + i..], lda)
                if (i < n - 1) {
                    (e[i], tau[i]) = dlarfg(n - i - 1, a[(i + 1) * lda + i], a[min(i + 2, n - 1) * lda + i..], lda)
                    a[(i + 1) * lda + i] = 1.0
                    nFloat64.dsymv(Lower, n - i - 1, 1.0, a[(i + 1) * lda + i + 1..], lda, a[(i + 1) * lda + i..], lda,
                        0.0, w[(i + 1) * ldw + i..], ldw)
                    nFloat64.dgemv(Trans, n - i - 1, i, 1.0, w[(i + 1) * ldw..], ldw, a[(i + 1) * lda + i..], lda, 0.0,
                        w[i..], ldw)
                    nFloat64.dgemv(NoTrans, n - i - 1, i, -1.0, a[(i + 1) * lda..], lda, w[i..], ldw, 1.0,
                        w[(i + 1) * ldw + i..], ldw)
                    nFloat64.dgemv(Trans, n - i - 1, i, 1.0, a[(i + 1) * lda..], lda, a[(i + 1) * lda + i..], lda, 0.0,
                        w[i..], ldw)
                    nFloat64.dgemv(NoTrans, n - i - 1, i, -1.0, w[(i + 1) * ldw..], ldw, w[i..], ldw, 1.0,
                        w[(i + 1) * ldw + i..], ldw)
                    nFloat64.dscal(n - i - 1, tau[i], w[(i + 1) * ldw + i..], ldw)
                    let alpha = -0.5 * tau[i] * nFloat64.ddot(n - i - 1, w[(i + 1) * ldw + i..], ldw,
                        a[(i + 1) * lda + i..], lda)
                    nFloat64.daxpy(n - i - 1, alpha, a[(i + 1) * lda + i..], lda, w[(i + 1) * ldw + i..], ldw)
                }
            }
        }
    }
}
