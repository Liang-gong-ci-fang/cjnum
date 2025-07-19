package cjnum.lapack.cjnum

import cjnum.blas.blas64.*
import cjnum.blas.*
import cjnum.lapack.*
import std.math.*

extend LFloat64Implementation {
    public func dsytrd(uplo: Uplo, n: Int64, a: Array<Float64>, lda: Int64, d: Array<Float64>, e: Array<Float64>,
        tau: Array<Float64>, work: Array<Float64>, lwork: Int64): Unit {
        if (uplo != Upper && uplo != Lower) {
            throw Exception(badUplo)
        } else if (n < 0) {
            throw Exception(nLT0)
        } else if (lda < max(1, n)) {
            throw Exception(badLdA)
        } else if (lwork < 1 && lwork != -1) {
            throw Exception(badLWork)
        } else if (work.size < max(1, lwork)) {
            throw Exception(shortWork)
        }

        if (n == 0) {
            work[0] = 1.0
            return
        }

        var nb = ilaenv(1, "DSYTRD", uplo.toString(), n, -1, -1, -1)
        let lworkopt = n * nb
        if (lwork == -1) {
            work[0] = Float64(lworkopt)
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

        var nx = n
        var iws = 1
        var ldwork: Int64
        if (1 < nb && nb < n) {
            nx = max(nb, ilaenv(3, "DSYTRD", uplo.toString(), n, -1, -1, -1))
            if (nx < n) {
                ldwork = nb
                iws = n * ldwork
                if (lwork < iws) {
                    nb = max(lwork / n, 1)
                    let nbmin = ilaenv(2, "DSYTRD", uplo.toString(), n, -1, -1, -1)
                    if (nb < nbmin) {
                        nx = n
                    }
                }
            } else {
                nx = n
            }
        } else {
            nb = 1
        }
        ldwork = nb

        if (uplo == Upper) {
            var i: Int64 = n - nb
            let kk = n - ((n - nx + nb - 1) / nb) * nb
            while (i >= kk) {
                dlatrd(uplo, i + nb, nb, a, lda, e, tau, work, ldwork)

                nFloat64.dsyr2k(uplo, NoTrans, i, nb, -1.0, a[i..], lda, work, ldwork, 1.0, a, lda)

                for (j in i..i + nb) {
                    a[(j - 1) * lda + j] = e[j - 1]
                    d[j] = a[j * lda + j]
                }
                i -= nb
            }
            dsytd2(uplo, kk, a, lda, d, e, tau)
        } else {
            var i: Int64 = 0
            while (i < n - nx) {
                dlatrd(uplo, n - i, nb, a[i * lda + i..], lda, e[i..], tau[i..], work, ldwork)

                nFloat64.dsyr2k(uplo, NoTrans, n - i - nb, nb, -1.0, a[(i + nb) * lda + i..], lda, work[nb * ldwork..],
                    ldwork, 1.0, a[(i + nb) * lda + i + nb..], lda)

                for (j in i..i + nb) {
                    a[(j + 1) * lda + j] = e[j]
                    d[j] = a[j * lda + j]
                }
                i += nb
            }
            dsytd2(uplo, n - i, a[i * lda + i..], lda, d[i..], e[i..], tau[i..])
        }
        work[0] = Float64(iws)
    }
}
