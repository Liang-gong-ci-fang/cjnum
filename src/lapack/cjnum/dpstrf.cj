package cjnum.lapack.cjnum

import cjnum.blas.blas64.*
import cjnum.blas.*
import cjnum.lapack.*
import std.math.*

extend LFloat64Implementation {
    public func dpstrf(uplo: Uplo, n: Int64, a: Array<Float64>, lda: Int64, piv: Array<Int64>, tol: Float64,
        work: Array<Float64>): (Int64, Bool) {
        if (uplo != Upper && uplo != Lower) {
            throw Exception(badUplo)
        }
        if (n < 0) {
            throw Exception(nLT0)
        }
        if (lda < max(1, n)) {
            throw Exception(badLdA)
        }

        if (n == 0) {
            return (0, true)
        }

        if (a.size < (n - 1) * lda + n) {
            throw Exception(shortA)
        }
        if (piv.size != n) {
            throw Exception(badLenPiv)
        }
        if (work.size < 2 * n) {
            throw Exception(shortWork)
        }

        let nb = ilaenv(1, "DPOTRF", uplo.toString(), n, -1, -1, -1)
        if (nb <= 1 || n <= nb) {
            return dpstf2(uplo, n, a, lda, piv, tol, work)
        }

        for (i in 0..piv[..n].size) {
            piv[i] = i
        }

        var pvt = 0
        var ajj = a[0]
        for (i in 1..n) {
            let aii = a[i * lda + i]
            if (aii > ajj) {
                pvt = i
                ajj = aii
            }
        }
        if (ajj <= 0.0 || ajj.isNaN()) {
            return (0, false)
        }

        var dstop = tol
        if (dstop < 0.0) {
            dstop = Float64(n) * dlamchE * ajj
        }

        let nFloat64 = nFloat64Implementation()
        let dots = work[..n]
        let work2 = work[n..2 * n]
        if (uplo == Upper) {
            for (k in 0..n : nb) {
                let jb = min(nb, n - k)
                for (i in k..n) {
                    dots[i] = 0.0
                }
                for (j in k..k + jb) {
                    for (i in j..n) {
                        if (j > k) {
                            let tmp = a[(j - 1) * lda + i]
                            dots[i] += tmp * tmp
                        }
                        work2[i] = a[i * lda + i] - dots[i]
                    }
                    if (j > 0) {
                        pvt = j
                        ajj = work2[pvt]
                        for (l in j + 1..n) {
                            let wl = work2[l]
                            if (wl > ajj) {
                                pvt = l
                                ajj = wl
                            }
                        }
                        if (ajj <= dstop || ajj.isNaN()) {
                            a[j * lda + j] = ajj
                            return (j, false)
                        }
                    }
                    if (j != pvt) {
                        a[pvt * lda + pvt] = a[j * lda + j]
                        nFloat64.dswap(j, a[j..], lda, a[pvt..], lda)
                        if (pvt < n - 1) {
                            nFloat64.dswap(n - pvt - 1, a[j * lda + (pvt + 1)..], 1, a[pvt * lda + (pvt + 1)..], 1)
                        }
                        nFloat64.dswap(pvt - j - 1, a[j * lda + (j + 1)..], 1, a[(j + 1) * lda + pvt..], lda)
                        (dots[j], dots[pvt]) = (dots[pvt], dots[j])
                        (piv[j], piv[pvt]) = (piv[pvt], piv[j])
                    }
                    ajj = sqrt(ajj)
                    a[j * lda + j] = ajj
                    if (j < n - 1) {
                        nFloat64.dgemv(Trans, j - k, n - j - 1, -1.0, a[k * lda + j + 1..], lda, a[k * lda + j..], lda,
                            1.0, a[j * lda + j + 1..], 1)
                        nFloat64.dscal(n - j - 1, 1.0 / ajj, a[j * lda + j + 1..], 1)
                    }
                }
                if (k + jb < n) {
                    let j = k + jb
                    nFloat64.dsyrk(Upper, Trans, n - j, jb, -1.0, a[k * lda + j..], lda, 1.0, a[j * lda + j..], lda)
                }
            }
        } else {
            for (k in 0..n : nb) {
                let jb = min(nb, n - k)
                for (i in k..n) {
                    dots[i] = 0.0
                }
                for (j in k..k + jb) {
                    for (i in j..n) {
                        if (j > k) {
                            let tmp = a[i * lda + (j - 1)]
                            dots[i] += tmp * tmp
                        }
                        work2[i] = a[i * lda + i] - dots[i]
                    }
                    if (j > 0) {
                        pvt = j
                        ajj = work2[pvt]
                        for (l in j + 1..n) {
                            let wl = work2[l]
                            if (wl > ajj) {
                                pvt = l
                                ajj = wl
                            }
                        }
                        if (ajj <= dstop || ajj.isNaN()) {
                            a[j * lda + j] = ajj
                            return (j, false)
                        }
                    }
                    if (j != pvt) {
                        a[pvt * lda + pvt] = a[j * lda + j]
                        nFloat64.dswap(j, a[j * lda..], 1, a[pvt * lda..], 1)
                        if (pvt < n - 1) {
                            nFloat64.dswap(n - pvt - 1, a[(pvt + 1) * lda + j..], lda, a[(pvt + 1) * lda + pvt..], lda)
                        }
                        nFloat64.dswap(pvt - j - 1, a[(j + 1) * lda + j..], lda, a[pvt * lda + (j + 1)..], 1)
                        (dots[j], dots[pvt]) = (dots[pvt], dots[j])
                        (piv[j], piv[pvt]) = (piv[pvt], piv[j])
                    }
                    ajj = sqrt(ajj)
                    a[j * lda + j] = ajj
                    if (j < n - 1) {
                        nFloat64.dgemv(NoTrans, n - j - 1, j - k, -1.0, a[(j + 1) * lda + k..], lda, a[j * lda + k..], 1,
                            1.0, a[(j + 1) * lda + j..], lda)
                        nFloat64.dscal(n - j - 1, 1.0 / ajj, a[(j + 1) * lda + j..], lda)
                    }
                }
                if (k + jb < n) {
                    let j = k + jb
                    nFloat64.dsyrk(Lower, NoTrans, n - j, jb, -1.0, a[j * lda + k..], lda, 1.0, a[j * lda + j..], lda)
                }
            }
        }
        return (n, true)
    }
}
