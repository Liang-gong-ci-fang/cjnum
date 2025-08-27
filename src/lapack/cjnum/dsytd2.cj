package cjnum.lapack.cjnum

import cjnum.blas.blas64.*
import cjnum.blas.*

extend LFloat64Implementation {
    public func dsytd2(uplo: Uplo, n: Int64, a: Array<Float64>, lda: Int64, d: Array<Float64>, e: Array<Float64>,
        tau: Array<Float64>): Unit {
        if (uplo != Upper && uplo != Lower) {
            throw Exception(badUplo)
        } else if (n < 0) {
            throw Exception(nLT0)
        } else if (lda < max(1, n)) {
            throw Exception(badLdA)
        }

        if (n == 0) {
            return
        }

        if (a.size < (n - 1) * lda + n) {
            throw Exception(shortA)
        } else if (d.size < n) {
            throw Exception(shortD)
        } else if (e.size < n - 1) {
            throw Exception(shortE)
        } else if (tau.size < n - 1) {
            throw Exception(shortTau)
        }

        let nFloat64 = nFloat64Implementation()

        if (uplo == Upper) {
            for (i in n - 2..=0 : -1) {
                var taui: Float64
                (a[i * lda + i + 1], taui) = dlarfg(i + 1, a[i * lda + i + 1], a[i + 1..], lda)
                e[i] = a[i * lda + i + 1]
                if (taui != 0.0) {
                    a[i * lda + i + 1] = 1.0

                    nFloat64.dsymv(uplo, i + 1, taui, a, lda, a[i + 1..], lda, 0.0, tau, 1)

                    let alpha = -0.5 * taui * nFloat64.ddot(i + 1, tau, 1, a[i + 1..], lda)
                    nFloat64.daxpy(i + 1, alpha, a[i + 1..], lda, tau, 1)

                    nFloat64.dsyr2(uplo, i + 1, -1.0, a[i + 1..], lda, tau, 1, a, lda)
                    a[i * lda + i + 1] = e[i]
                }
                d[i + 1] = a[(i + 1) * lda + i + 1]
                tau[i] = taui
            }
            d[0] = a[0]
            return
        }
        for (i in 0..n - 1) {
            var taui: Float64
            (a[(i + 1) * lda + i], taui) = dlarfg(n - i - 1, a[(i + 1) * lda + i], a[min(i + 2, n - 1) * lda + i..], lda
            )
            e[i] = a[(i + 1) * lda + i]
            if (taui != 0.0) {
                a[(i + 1) * lda + i] = 1.0

                nFloat64.dsymv(uplo, n - i - 1, taui, a[(i + 1) * lda + i + 1..], lda, a[(i + 1) * lda + i..], lda, 0.0,
                    tau[i..], 1)

                let alpha = -0.5 * taui * nFloat64.ddot(n - i - 1, tau[i..], 1, a[(i + 1) * lda + i..], lda)
                nFloat64.daxpy(n - i - 1, alpha, a[(i + 1) * lda + i..], lda, tau[i..], 1)

                nFloat64.dsyr2(uplo, n - i - 1, -1.0, a[(i + 1) * lda + i..], lda, tau[i..], 1, a[(i + 1) * lda + i + 1..],
                    lda)
                a[(i + 1) * lda + i] = e[i]
            }
            d[i] = a[i * lda + i]
            tau[i] = taui
        }
        d[n - 1] = a[(n - 1) * lda + n - 1]
    }
}
