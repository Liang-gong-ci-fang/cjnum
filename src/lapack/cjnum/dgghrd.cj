package cjnum.lapack.cjnum

import cjnum.blas.blas64.*
import cjnum.blas.*
import cjnum.lapack.*
import std.math.*

extend LFloat64Implementation {
    public func dgghrd(compq: OrthoComp, compz: OrthoComp, n: Int, ilo: Int, ihi: Int, a: Array<Float64>, lda: Int,
        b: Array<Float64>, ldb: Int, q: Array<Float64>, ldq: Int, z: Array<Float64>, ldz: Int): Unit {
        if (compq != OrthoNone && compq != OrthoExplicit && compq != OrthoPostmul) {
            throw Exception(badOrthoComp)
        } else if (compz != OrthoNone && compz != OrthoExplicit && compz != OrthoPostmul) {
            throw Exception(badOrthoComp)
        } else if (n < 0) {
            throw Exception(nLT0)
        } else if (ilo < 0 || max(0, n - 1) < ilo) {
            throw Exception(badIlo)
        } else if (ihi < min(ilo, n - 1) || n <= ihi) {
            throw Exception(badIhi)
        } else if (lda < max(1, n)) {
            throw Exception(badLdA)
        } else if (ldb < max(1, n)) {
            throw Exception(badLdB)
        } else if ((compq != OrthoNone && ldq < n) || ldq < 1) {
            throw Exception(badLdQ)
        } else if ((compz != OrthoNone && ldz < n) || ldz < 1) {
            throw Exception(badLdZ)
        }

        if (n == 0) {
            return
        }

        if (a.size < (n - 1) * lda + n) {
            throw Exception(shortA)
        } else if (b.size < (n - 1) * ldb + n) {
            throw Exception(shortB)
        } else if (compq != OrthoNone && q.size < (n - 1) * ldq + n) {
            throw Exception(shortQ)
        } else if (compz != OrthoNone && z.size < (n - 1) * ldz + n) {
            throw Exception(shortZ)
        }

        if (compq == OrthoExplicit) {
            dlaset(All, n, n, 0.0, 1.0, q, ldq)
        }
        if (compz == OrthoExplicit) {
            dlaset(All, n, n, 0.0, 1.0, z, ldz)
        }

        if (n == 1) {
            return
        }

        for (i in 1..n) {
            for (j in 0..i) {
                b[i * ldb + j] = 0.0
            }
        }
        let nFloat64 = nFloat64Implementation()

        for (jcol in ilo..=ihi - 2) {
            for (jrow in ihi..=jcol + 2 : -1) {
                var (c, s) = (0.0, 0.0)
                (c, s, a[(jrow - 1) * lda + jcol]) = dlartg(a[(jrow - 1) * lda + jcol], a[jrow * lda + jcol])
                a[jrow * lda + jcol] = 0.0

                nFloat64.drot(n - jcol - 1, a[(jrow - 1) * lda + jcol + 1..], 1, a[jrow * lda + jcol + 1..], 1, c, s)
                nFloat64.drot(n + 2 - jrow - 1, b[(jrow - 1) * ldb + jrow - 1..], 1, b[jrow * ldb + jrow - 1..], 1, c, s
                )

                if (compq != OrthoNone) {
                    nFloat64.drot(n, q[jrow - 1..], ldq, q[jrow..], ldq, c, s)
                }

                (c, s, b[jrow * ldb + jrow]) = dlartg(b[jrow * ldb + jrow], b[jrow * ldb + jrow - 1])
                b[jrow * ldb + jrow - 1] = 0.0

                nFloat64.drot(ihi + 1, a[jrow..], lda, a[jrow - 1..], lda, c, s)
                nFloat64.drot(jrow, b[jrow..], ldb, b[jrow - 1..], ldb, c, s)

                if (compz != OrthoNone) {
                    nFloat64.drot(n, z[jrow..], ldz, z[jrow - 1..], ldz, c, s)
                }
            }
        }
    }
}
