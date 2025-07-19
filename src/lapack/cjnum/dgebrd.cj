package cjnum.lapack.cjnum

import cjnum.blas.blas64.*
import cjnum.blas.*
import std.math.*

extend LFloat64Implementation {
    public func dgebrd(m: Int64, n: Int64, a: Array<Float64>, lda: Int64, d: Array<Float64>, e: Array<Float64>,
        tauQ: Array<Float64>, tauP: Array<Float64>, work: Array<Float64>, lwork: Int64): Unit {
        if (m < 0) {
            throw Exception(mLT0)
        }
        if (n < 0) {
            throw Exception(nLT0)
        }
        if (lda < max(1, n)) {
            throw Exception(badLdA)
        }
        if (lwork < max(1, max(m, n)) && lwork != -1) {
            throw Exception(badLWork)
        }
        if (work.size < max(1, lwork)) {
            throw Exception(shortWork)
        }

        let minmn = min(m, n)
        if (minmn == 0) {
            work[0] = 1.0
            return
        }

        var nb = ilaenv(1, "DGEBRD", " ", m, n, -1, -1)
        let lwkopt = (m + n) * nb
        if (lwork == -1) {
            work[0] = Float64(lwkopt)
            return
        }

        if (a.size < (m - 1) * lda + n) {
            throw Exception(shortA)
        }
        if (d.size < minmn) {
            throw Exception(shortD)
        }
        if (e.size < minmn - 1) {
            throw Exception(shortE)
        }
        if (tauQ.size < minmn) {
            throw Exception(shortTauQ)
        }
        if (tauP.size < minmn) {
            throw Exception(shortTauP)
        }

        var nx = minmn
        var ws = max(m, n)
        if (1 < nb && nb < minmn) {
            nx = max(nb, ilaenv(3, "DGEBRD", " ", m, n, -1, -1))
            if (nx < minmn) {
                let ws = (m + n) * nb
                if (lwork < ws) {
                    let nbmin = ilaenv(2, "DGEBRD", " ", m, n, -1, -1)
                    if (lwork >= (m + n) * nbmin) {
                        nb = lwork / (m + n)
                    } else {
                        nb = minmn
                        nx = minmn
                    }
                }
            }
        }
        let nFloat64 = nFloat64Implementation()
        let ldworkx = nb
        let ldworky = nb
        var i: Int64 = 0
        while (i < minmn - nx) {
            let x = work[..m * ldworkx]
            let y = work[m * ldworkx..]
            dlabrd(m - i, n - i, nb, a[i * lda + i..], lda, d[i..], e[i..], tauQ[i..], tauP[i..], x, ldworkx, y, ldworky
            )

            nFloat64.dgemm(NoTrans, Trans, m - i - nb, n - i - nb, nb, -1.0, a[(i + nb) * lda + i..], lda,
                y[nb * ldworky..], ldworky, 1.0, a[(i + nb) * lda + i + nb..], lda)

            nFloat64.dgemm(NoTrans, NoTrans, m - i - nb, n - i - nb, nb, -1.0, x[nb * ldworkx..], ldworkx,
                a[i * lda + i + nb..], lda, 1.0, a[(i + nb) * lda + i + nb..], lda)

            if (m >= n) {
                for (j in i..i + nb) {
                    a[j * lda + j] = d[j]
                    a[j * lda + j + 1] = e[j]
                }
            } else {
                for (j in i..i + nb) {
                    a[j * lda + j] = d[j]
                    a[(j + 1) * lda + j] = e[j]
                }
            }
            i += nb
        }
        dgebd2(m - i, n - i, a[i * lda + i..], lda, d[i..], e[i..], tauQ[i..], tauP[i..], work)
        work[0] = Float64(ws)
    }
}
