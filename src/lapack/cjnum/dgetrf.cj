package cjnum.lapack.cjnum

import std.math.*
import cjnum.lapack.*
import cjnum.blas.*
import cjnum.blas.blas64.*

extend LFloat64Implementation {
    public func dgetrf(m: Int64, n: Int64, a: Array<Float64>, lda: Int64, ipiv: Array<Int64>): Bool {
        // TODO
        throw Exception("not implemented")
        // let mn = min(m, n)
        // if (m < 0) {
        //     throw Exception(mLT0)
        // }
        // if (n < 0) {
        //     throw Exception(nLT0)
        // }
        // if (lda < max(1, n)) {
        //     throw Exception(badLdA)
        // }

        // if (mn == 0) {
        //     return true
        // }

        // if (a.size < (m - 1) * lda + n) {
        //     throw Exception(shortA)
        // }
        // if (ipiv.size < mn) {
        //     throw Exception(badLenIpiv)
        // }

        // let nFloat64 = nFloat64Implementation()

        // let nb = ilaenv(1, "DGETRF", " ", m, n, -1, -1)
        // if (nb <= 1 || mn <= nb) {
        //     return dgetf2(m, n, a, lda, ipiv)
        // }
        // var ok = true
        // var j = 0
        // while (j < mn) {
        //     let jb = min(mn - j, nb)
        //     let blockOk = dgetf2(m - j, jb, a[j * lda + j..], lda, ipiv[j..j + jb])
        //     if (!blockOk) {
        //         ok = false
        //     }
        //     for (i in j..min(m - 1, j + jb - 1)) {
        //         ipiv[i] = j + ipiv[i]
        //     }
        //     dlaswp(j, a, lda, j, j + jb - 1, ipiv[..j + jb], 1)
        //     if (j + jb < n) {
        //         dlaswp(n - j - jb, a[j + jb..], lda, j, j + jb - 1, ipiv[..j + jb], 1)
        //         nFloat64.dtrsm(Left, Lower, NoTrans, UnitDiag,
		// 		jb, n - j - jb, 1.0,
		// 		a[j * lda + j..], lda,
		// 		a[j * lda + j + jb..], lda)
        //         if (j + jb < m ) {
        //             nFloat64.dgemm(NoTrans, NoTrans, m - j - jb, n - j - jb, jb, -1.0,
		// 			a[(j + jb) * lda + j..], lda,
		// 			a[j * lda + j + jb..], lda,
		// 			1.0, a[(j + jb) * lda + j + jb..], lda)
        //         }
        //     }
        //     j += nb
        // }
        // return ok
    }
}