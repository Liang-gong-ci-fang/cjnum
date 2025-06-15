package cjnum.lapack.cjnum

import cjnum.blas.blas64.*
import cjnum.blas.*
import std.math.*

extend LFloat64Implementation {
    public func dlahr2(n: Int64, k: Int64, nb: Int64, a: Array<Float64>, lda: Int64, tau: Array<Float64>,
        t: Array<Float64>, ldt: Int64, y: Array<Float64>, ldy: Int64): Unit {
        if (n < 0) {
            throw Exception(nLT0)
        } else if (k < 0) {
            throw Exception(kLT0)
        } else if (nb < 0) {
            throw Exception(nbLT0)
        } else if (nb > n) {
            throw Exception(nbGTN)
        } else if (lda < max(1, n - k + 1)) {
            throw Exception(badLdA)
        } else if (ldt < max(1, nb)) {
            throw Exception(badLdT)
        } else if (ldy < max(1, nb)) {
            throw Exception(badLdY)
        }

        if (n < 0) {
            return
        }

        if (a.size < (n - 1) * lda + n - k + 1) {
            throw Exception(shortA)
        } else if (tau.size < nb) {
            throw Exception(shortTau)
        } else if (t.size < (nb - 1) * ldt + nb) {
            throw Exception(shortT)
        } else if (y.size < (n - 1) * ldy + nb) {
            throw Exception(shortY)
        }

        if (n == 1) {
            return
        }

        let nFloat64 = nFloat64Implementation()
        var ei: Float64 = 0.0
        for (i in 0..nb) {
            if (i > 0) {
                nFloat64.dgemv(NoTrans, n - k, i, -1.0, y[k * ldy..], ldy, a[(k + i - 1) * lda..], 1, 1.0,
                    a[k * lda + i..], lda)

                nFloat64.dcopy(i, a[k * lda + i..], lda, t[nb - 1..], ldt)
                nFloat64.dtrmv(Lower, Trans, UnitDiag, i, a[k * lda..], lda, t[nb - 1..], ldt)

                nFloat64.dgemv(Trans, n - k - i, i, 1.0, a[(k + i) * lda..], lda, a[(k + i) * lda + i..], lda, 1.0,
                    t[nb - 1..], ldt)

                nFloat64.dtrmv(Upper, Trans, NonUnitDiag, i, t, ldt, t[nb - 1..], ldt)

                nFloat64.dgemv(NoTrans, n - k - i, i, -1.0, a[(k + i) * lda..], lda, t[nb - 1..], ldt, 1.0,
                    a[(k + i) * lda + i..], lda)

                nFloat64.dtrmv(Lower, NoTrans, UnitDiag, i, a[k * lda..], lda, t[nb - 1..], ldt)
                nFloat64.daxpy(i, -1.0, t[nb - 1..], ldt, a[k * lda + i..], lda)

                a[(k + i - 1) * lda + i - 1] = ei
            }

            (ei, tau[i]) = dlarfg(n - k - i, a[(k + i) * lda + i], a[min(k + i + 1, n - 1) * lda + i..], lda)
            a[(k + i) * lda + i] = 1.0

            nFloat64.dgemv(NoTrans, n - k, n - k - i, 1.0, a[k * lda + i + 1..], lda, a[(k + i) * lda + i..], lda, 0.0,
                y[k * ldy + i..], ldy)
            nFloat64.dgemv(Trans, n - k - i, i, 1.0, a[(k + i) * lda..], lda, a[(k + i) * lda + i..], lda, 0.0, t[i..],
                ldt)
            nFloat64.dgemv(NoTrans, n - k, i, -1.0, y[k * ldy..], ldy, t[i..], ldt, 1.0, y[k * ldy + i..], ldy)
            nFloat64.dscal(n - k, tau[i], y[k * ldy + i..], ldy)

            nFloat64.dscal(i, -tau[i], t[i..], ldt)
            nFloat64.dtrmv(Upper, NoTrans, NonUnitDiag, i, t, ldt, t[i..], ldt)

            t[i * ldt + i] = tau[i]
        }
        a[(k + nb - 1) * lda + nb - 1] = ei

        dlacpy(All, k, nb, a[1..], lda, y, ldy)
        nFloat64.dtrmm(Right, Lower, NoTrans, UnitDiag, k, nb, 1.0, a[k * lda..], lda, y, ldy)
        if (n > k + nb) {
            nFloat64.dgemm(NoTrans, NoTrans, k, nb, n - k - nb, 1.0, a[1 + nb..], lda, a[(k + nb) * lda..], lda, 1.0, y,
                ldy)
        }
        nFloat64.dtrmm(Right, Upper, NoTrans, NonUnitDiag, k, nb, 1.0, t, ldt, y, ldy)
    }
}
