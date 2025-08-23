package cjnum.lapack.cjnum

import std.math.*
import cjnum.blas.*
import cjnum.blas.blas64.*
import cjnum.lapack.*

extend LFloat64Implementation {
    public func dgetri(n: Int64, a: Array<Float64>, lda: Int64, ipiv: Array<Int64>, work: Array<Float64>, lwork: Int64): Bool {
    // public func dgetri(_: Int64, _: Array<Float64>, _: Int64, _: Array<Int64>, _: Array<Float64>, _: Int64): Bool {
    //     // TODO
    //     throw Exception("not implemented")
        let iws = max(1, n)
        if (n < 0) {
            throw Exception(nLT0)
        }
        if (lda < max(1, n)) {
            throw Exception(badLdA)
        }
        if (lwork < iws && lwork != -1) {
            throw Exception(badLWork)
        }
        if (work.size < max(1, lwork)) {
            throw Exception(shortWork)
        }

        if (n == 0) {
            work[0] = 1.0
            return true
        }

        var nb = ilaenv(1, "DGETRI", " ", n, -1, -1, -1)
        if (lwork == -1) {
            work[0] = Float64(n * nb)
            return true
        }

        if (a.size < (n - 1) * lda + n) {
            throw Exception(shortA)
        }
        if (ipiv.size != n) {
            throw Exception(badLenIpiv)
        }

        let ok = dtrtri(Upper, NonUnitDiag, n, a, lda)
        if (!ok) {
            return false
        }

        var nbmin = 2
        if (1 < nb && nb < n) {
            let iws = max(n * nb, 1)
            if (lwork < iws) {
                nb = lwork / n
                nbmin = max(2, ilaenv(2, "DGETRI", " ", n, -1, -1, -1))
            }
        }
        let ldwork = nb

        let nFloat64 = nFloat64Implementation()

        if (nb < nbmin || n <= nb) {
            var j = n - 1
            while (j >= 0) {
                for (i in j + 1..n) {
                    work[i] = a[i * lda + j]
                    a[i * lda + j] = 0.0
                }
                if (j < n - 1) {
                    nFloat64.dgemv(NoTrans, n, n - j - 1, -1.0, a[(j + 1)..], lda, work[(j + 1)..], 1, 1.0, a[j..], lda)
                }
                j--
            }
        } else {
            let nn = ((n - 1) / nb) * nb
            var j = nn
            while (j >= 0) {
                let jb = min(nb, n - j)
                for (jj in j ..j + jb) {
                    for (i in jj + 1..n) {
                        work[i * ldwork + (jj - j)] = a[i * lda + jj]
                        a[i * lda + jj] = 0.0
                    }
                }
                if (j + jb < n) {
                    nFloat64.dgemm(NoTrans, NoTrans, n, jb, n - j - jb, -1.0, a[(j + jb)..], lda, work[(j + jb) * ldwork..], ldwork, 1.0, a[j..], lda)
                }
                nFloat64.dtrsm(Right, Lower, NoTrans, UnitDiag, n, jb, 1.0, work[j * ldwork..], ldwork, a[j..], lda)
                j -= nb
            }
        }

        var j = n - 2
        while (j >= 0) {
            let jp = ipiv[j]
            if (jp != j) {
                nFloat64.dswap(n, a[j..], lda, a[jp..], lda)
            }
            j--
        }
        work[0] = Float64(iws)

        return true
    }
}